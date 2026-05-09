import gmsh
import numpy as np
import tqdm.autonotebook
import petsc4py
import time

petsc4py.init()
from petsc4py import PETSc

from mpi4py import MPI
from basix.ufl import element
from dolfinx.fem import (
    Constant, Function, functionspace,
    dirichletbc, extract_function_spaces,
    form, locate_dofs_topological,
)
from dolfinx.fem.petsc import (
    apply_lifting, assemble_matrix, assemble_vector,
    create_vector, set_bc,
)
from dolfinx.io import VTXWriter, gmsh as gmshio
from ufl import (
    FacetNormal, TestFunction, TrialFunction,
    div, dot, dx, ds, inner, lhs, nabla_grad, rhs, sym, Identity,
)

from pathlib import Path

BASE_DIR = Path(__file__).parent

mesh_comm = MPI.COMM_WORLD
model_rank = 0

inlet_marker, outlet_marker, wall_marker, car_marker, air_marker = 1, 2, 3, 4, 5
gmsh.initialize()
gmsh.merge(str(BASE_DIR / "pickup_mock.msh"))
mesh_data = gmshio.model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=3)
mesh = mesh_data.mesh
ft = mesh_data.facet_tags
ct = mesh_data.cell_tags
gmsh.finalize()

t = 0.0
T = 3.0
dt = 0.01
num_steps = int(T / dt)
U_inf = 16.0
k = Constant(mesh, dt)
mu = Constant(mesh, 1.8e-5)
rho = Constant(mesh, 1.225)

v_cg2 = element("Lagrange", mesh.basix_cell(), 2, shape=(mesh.geometry.dim,))
s_cg1 = element("Lagrange", mesh.basix_cell(), 1)
V = functionspace(mesh, v_cg2)
Q = functionspace(mesh, s_cg1)


def inlet_velocity(x):
    values = np.zeros((3, x.shape[1]), dtype=float)
    values[0] = U_inf
    return values


# Граничные условия
u_inlet = Function(V)
u_inlet.interpolate(inlet_velocity)
bcu_inflow = dirichletbc(u_inlet, locate_dofs_topological(V, 2, ft.find(inlet_marker)))

u_nonslip = np.array((0,) * mesh.geometry.dim, dtype=float)
bcu_walls = dirichletbc(u_nonslip, locate_dofs_topological(V, 2, ft.find(wall_marker)), V)
bcu_car = dirichletbc(u_nonslip, locate_dofs_topological(V, 2, ft.find(car_marker)), V)

bcu = [bcu_inflow, bcu_walls, bcu_car]

bcp_outlet = dirichletbc(0.0, locate_dofs_topological(Q, 2, ft.find(outlet_marker)), Q)

bcp = [bcp_outlet]

# Уравнения
u = TrialFunction(V)
v = TestFunction(V)
u_ = Function(V)
u_s = Function(V)
u_n = Function(V)
p = TrialFunction(Q)
q = TestFunction(Q)
p_ = Function(Q)
p_n = Function(Q)
n = FacetNormal(mesh)
f = Constant(mesh, np.zeros(3))


def epsilon(u):
    return sym(nabla_grad(u))


def sigma(u, p):
    return 2 * mu * epsilon(u) - p * Identity(3)


F1 = rho * dot((u - u_n) / k, v) * dx
F1 += rho * dot(dot(u_n, nabla_grad(u_n)), v) * dx
F1 += inner(sigma((u + u_n) / 2, p_n), epsilon(v)) * dx
F1 += dot(p_n * n, v) * ds
F1 -= dot(mu * dot(nabla_grad((u + u_n) / 2), n), v) * ds
F1 -= dot(f, v) * dx
a1 = form(lhs(F1))
L1 = form(rhs(F1))
A1 = assemble_matrix(a1, bcs=bcu)
A1.assemble()
b1 = create_vector(extract_function_spaces(L1))

a2 = form(dot(nabla_grad(p), nabla_grad(q)) * dx)
L2 = form(dot(nabla_grad(p_n), nabla_grad(q)) * dx - (rho / k) * div(u_s) * q * dx)
A2 = assemble_matrix(a2, bcs=bcp)
A2.assemble()
b2 = create_vector(extract_function_spaces(L2))

a3 = form(dot(u, v) * dx)
L3 = form(dot(u_s, v) * dx - k * dot(nabla_grad(p_ - p_n), v) * dx)
A3 = assemble_matrix(a3)
A3.assemble()
b3 = create_vector(extract_function_spaces(L3))

# Солверы
solver1 = PETSc.KSP().create(mesh.comm)
solver1.setOperators(A1)
solver1.setType("bcgs")
pc1 = solver1.getPC()
pc1.setType("ilu")
solver1.setTolerances(rtol=1e-6, atol=1e-10)

solver2 = PETSc.KSP().create(mesh.comm)
solver2.setOperators(A2)
solver2.setType("cg")
pc2 = solver2.getPC()
try:
    pc2.setType("hypre")
    pc2.setHYPREType("boomeramg")
except:
    pc2.setType("gamg")
solver2.setTolerances(rtol=1e-8, atol=1e-12)

solver3 = PETSc.KSP().create(mesh.comm)
solver3.setOperators(A3)
solver3.setType("cg")
pc3 = solver3.getPC()
pc3.setType("sor")
solver3.setTolerances(rtol=1e-6, atol=1e-10)

# Расчёт
folder = Path(str(BASE_DIR / "results"))
folder.mkdir(exist_ok=True, parents=True)
vtx_u = VTXWriter(mesh.comm, folder / "u.bp", [u_], engine="BP4")
vtx_p = VTXWriter(mesh.comm, folder / "p.bp", [p_], engine="BP4")
vtx_u.write(t)
vtx_p.write(t)
progress = tqdm.autonotebook.tqdm(desc="Solving PDE", total=num_steps)
for i in range(num_steps):
    progress.update(1)
    t += dt

    with b1.localForm() as loc:
        loc.set(0)
    assemble_vector(b1, L1)
    apply_lifting(b1, [a1], [bcu])
    b1.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b1, bcu)
    solver1.solve(b1, u_s.x.petsc_vec)
    u_s.x.scatter_forward()

    with b2.localForm() as loc:
        loc.set(0)
    assemble_vector(b2, L2)
    apply_lifting(b2, [a2], [bcp])
    b2.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b2, bcp)
    solver2.solve(b2, p_.x.petsc_vec)
    p_.x.scatter_forward()

    with b3.localForm() as loc:
        loc.set(0)
    assemble_vector(b3, L3)
    b3.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    solver3.solve(b3, u_.x.petsc_vec)
    u_.x.scatter_forward()

    vtx_u.write(t)
    vtx_p.write(t)

    with (
        u_.x.petsc_vec.localForm() as loc_,
        u_n.x.petsc_vec.localForm() as loc_n,
    ):
        loc_.copy(loc_n)
    with (
        p_.x.petsc_vec.localForm() as loc_p,
        p_n.x.petsc_vec.localForm() as loc_pn
    ):
        loc_p.copy(loc_pn)
