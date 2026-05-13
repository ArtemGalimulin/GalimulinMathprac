import gmsh
import numpy as np
import tqdm.autonotebook
import petsc4py

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

# mpirun -n 4 python3 navier-stokes.py

BASE_DIR = Path(__file__).parent

mesh_comm = MPI.COMM_WORLD
model_rank = 0

inlet_marker, outlet_marker, wall_marker, floor_marker, car_marker, air_marker = 1, 2, 3, 4, 5, 6
gmsh.initialize()
gmsh.merge(str(BASE_DIR / "light_03_relative.msh"))
mesh_data = gmshio.model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=3)
mesh = mesh_data.mesh
ft = mesh_data.facet_tags
ct = mesh_data.cell_tags
gmsh.finalize()

t = 0.0
T = 3.0
t_acc = 0.1
dt = 0.005
num_steps = int(T / dt)
U_inf = 10.0
k = Constant(mesh, dt)
mu = Constant(mesh, 0.1)
rho = Constant(mesh, 1.0)

v_cg2 = element("Lagrange", mesh.basix_cell(), 2, shape=(mesh.geometry.dim,))
s_cg1 = element("Lagrange", mesh.basix_cell(), 1)
V = functionspace(mesh, v_cg2)
Q = functionspace(mesh, s_cg1)


class InletVelocity:
    def __init__(self):
        self.t = 0.0

    def __call__(self, x):
        ramp = min(self.t / t_acc, 1.0)
        values = np.zeros((3, x.shape[1]), dtype=float)
        values[1] = U_inf * ramp
        return values


# Граничные условия
inlet_velocity = InletVelocity()
u_inlet = Function(V)
u_inlet.interpolate(inlet_velocity)
bcu_inflow = dirichletbc(u_inlet, locate_dofs_topological(V, 2, ft.find(inlet_marker)))
bcu_floor = dirichletbc(u_inlet, locate_dofs_topological(V, 2, ft.find(floor_marker)))
bcu_walls = dirichletbc(u_inlet, locate_dofs_topological(V, 2, ft.find(wall_marker)))

u_nonslip = np.array((0,) * mesh.geometry.dim, dtype=float)
# bcu_floor = dirichletbc(u_nonslip, locate_dofs_topological(V, 2, ft.find(floor_marker)), V)
# bcu_walls = dirichletbc(u_nonslip, locate_dofs_topological(V, 2, ft.find(wall_marker)), V)
bcu_car = dirichletbc(u_nonslip, locate_dofs_topological(V, 2, ft.find(car_marker)), V)

bcu = [bcu_inflow, bcu_floor, bcu_walls, bcu_car]

bcp_outlet = dirichletbc(PETSc.ScalarType(0.0), locate_dofs_topological(Q, 2, ft.find(outlet_marker)), Q)

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
# вот здесь в статье плотность отсутствует, а она нужна
L3 = form(dot(u_s, v) * dx - (k / rho) * dot(nabla_grad(p_ - p_n), v) * dx)
A3 = assemble_matrix(a3, bcs=bcu)
A3.assemble()
b3 = create_vector(extract_function_spaces(L3))

# Солверы
solver1 = PETSc.KSP().create(mesh.comm)
solver1.setOperators(A1)
solver1.setType("bcgs")  # bicgstab
pc1 = solver1.getPC()
pc1.setType("hypre")
pc1.setHYPREType("boomeramg")
solver1.setTolerances(rtol=1e-6, atol=1e-10)

# Шаг 2: давление
solver2 = PETSc.KSP().create(mesh.comm)
solver2.setOperators(A2)
solver2.setType("bcgs")  # bicgstab
pc2 = solver2.getPC()
pc2.setType("hypre")
pc2.setHYPREType("boomeramg")
solver2.setTolerances(rtol=1e-8, atol=1e-12)

# Шаг 3: коррекция скорости (CG + SOR)
solver3 = PETSc.KSP().create(mesh.comm)
solver3.setOperators(A3)
solver3.setType("cg")  # conjugate gradient
pc3 = solver3.getPC()
pc3.setType("sor")  # successive over-relaxation
solver3.setTolerances(rtol=1e-8, atol=1e-10)

# Расчёт
folder = Path(str(BASE_DIR / "results"))
folder.mkdir(exist_ok=True, parents=True)
vtx_u = VTXWriter(mesh.comm, folder / "u5.bp", [u_], engine="BP4")
vtx_p = VTXWriter(mesh.comm, folder / "p5.bp", [p_], engine="BP4")
vtx_u.write(t)
vtx_p.write(t)
progress = tqdm.autonotebook.tqdm(desc="Solving PDE", total=num_steps)
for i in range(num_steps):
    progress.update(1)
    t += dt
    inlet_velocity.t = t
    u_inlet.interpolate(inlet_velocity)

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
    apply_lifting(b3, [a3], [bcu])
    b3.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b3, bcu)
    solver3.solve(b3, u_.x.petsc_vec)
    u_.x.scatter_forward()

    if i % 10 == 0:
        vtx_u.write(t)
        vtx_p.write(t)

        r1 = solver1.getResidualNorm()
        r2 = solver2.getResidualNorm()
        it1 = solver1.getIterationNumber()
        it2 = solver2.getIterationNumber()

        u_norm = u_.x.petsc_vec.norm()

        print(f"t={t:.3f} | u_norm={u_norm:.3e} | "
              f"KSP1: {it1}iter res={r1:.2e} | "
              f"KSP2: {it2}iter res={r2:.2e}")

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

    u_n.x.petsc_vec.ghostUpdate(addv=PETSc.InsertMode.INSERT_VALUES, mode=PETSc.ScatterMode.FORWARD)
    p_n.x.petsc_vec.ghostUpdate(addv=PETSc.InsertMode.INSERT_VALUES, mode=PETSc.ScatterMode.FORWARD)