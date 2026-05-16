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
    form, locate_dofs_topological, Expression,
)
from dolfinx.fem.petsc import (
    apply_lifting, assemble_matrix, assemble_vector,
    create_vector, set_bc
)
from dolfinx.fem import assemble_scalar
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
if mesh_comm.rank == 0:
    gmsh.model.add("schaefer_turek_3d")

    L, W, H = 2.5, 0.41, 0.41
    x_cyl, y_cyl, r_cyl = 0.5, 0.20, 0.05

    box = gmsh.model.occ.addBox(0, 0, 0, L, W, H)
    cyl = gmsh.model.occ.addCylinder(x_cyl, y_cyl, 0, 0, 0, H, r_cyl)
    fluid, _ = gmsh.model.occ.cut([(3, box)], [(3, cyl)])
    gmsh.model.occ.synchronize()

    # 2. Поиск и маркировка поверхностей через Центр Масс
    surfaces = gmsh.model.getEntities(dim=2)
    inlet_surfs, outlet_surfs, wall_surfs, floor_surfs, car_surfs = [], [], [], [], []
    for dim, tag in surfaces:
        com = gmsh.model.occ.getCenterOfMass(dim, tag)
        if np.isclose(com[0], 0.0):
            inlet_surfs.append(tag)
        elif np.isclose(com[0], L):
            outlet_surfs.append(tag)
        elif np.isclose(com[0], x_cyl) and np.isclose(com[1], y_cyl):
            car_surfs.append(tag)
        elif np.isclose(com[2], 0.0) or np.isclose(com[2], H):
            floor_surfs.append(tag)
        else:
            wall_surfs.append(tag)

    gmsh.model.addPhysicalGroup(2, inlet_surfs, inlet_marker)
    gmsh.model.addPhysicalGroup(2, outlet_surfs, outlet_marker)
    gmsh.model.addPhysicalGroup(2, wall_surfs, wall_marker)
    gmsh.model.addPhysicalGroup(2, floor_surfs, floor_marker)
    gmsh.model.addPhysicalGroup(2, car_surfs, car_marker)
    gmsh.model.addPhysicalGroup(3, [fluid[0][1]], air_marker)

    res_min = r_cyl / 3.5

    distance_field = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance_field, "FacesList", car_surfs)

    threshold_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
    gmsh.model.mesh.field.setNumber(threshold_field, "LcMin", res_min)
    gmsh.model.mesh.field.setNumber(threshold_field, "LcMax", 0.04)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", r_cyl)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", 3 * r_cyl)

    min_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", [threshold_field])
    gmsh.model.mesh.field.setAsBackgroundMesh(min_field)

    gmsh.option.setNumber("Mesh.Algorithm", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    gmsh.option.setNumber("Mesh.RecombineAll", 0)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(2)

    gmsh.model.mesh.optimize("HighOrder")
    gmsh.model.mesh.generate(3)

    gmsh.model.mesh.setOrder(2)

    gmsh.model.mesh.optimize("HighOrder")

mesh_data = gmshio.model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=3)
mesh = mesh_data.mesh
ft = mesh_data.facet_tags
ct = mesh_data.cell_tags
gmsh.finalize()

t = 0.0
T = 4.0
t_acc = 0.2
dt = 0.005
num_steps = int(T / dt)
save_freq = 10
U_m = 0.45
U_mean = 0.20
H_channel = 0.41
D_cylinder = 0.1
A_ref = D_cylinder * H_channel

k = Constant(mesh, dt)
mu = Constant(mesh, 0.001)
rho = Constant(mesh, 1.0)

v_cg2 = element("Lagrange", mesh.basix_cell(), 2, shape=(mesh.geometry.dim,))
s_cg1 = element("Lagrange", mesh.basix_cell(), 1)
V = functionspace(mesh, v_cg2)
Q = functionspace(mesh, s_cg1)


class InletVelocity:
    def __init__(self):
        self.t = 0.0

    def __call__(self, x):
        # ramp = min(self.t / t_acc, 1.0)
        ramp = 1
        values = np.zeros((3, x.shape[1]), dtype=float)
        values[0] = 16.0 * U_m * x[1] * x[2] * (H_channel - x[1]) * (H_channel - x[2]) / (H_channel**4) * ramp
        return values

# Граничные условия
inlet_velocity = InletVelocity()
u_inlet = Function(V)
u_inlet.interpolate(inlet_velocity)

u_nonslip = np.array((0,) * mesh.geometry.dim, dtype=float)

bcu_inflow = dirichletbc(u_inlet, locate_dofs_topological(V, 2, ft.find(inlet_marker)))
bcu_floor = dirichletbc(u_nonslip, locate_dofs_topological(V, 2, ft.find(floor_marker)), V)
bcu_walls = dirichletbc(u_nonslip, locate_dofs_topological(V, 2, ft.find(wall_marker)), V)
bcu_car = dirichletbc(u_nonslip, locate_dofs_topological(V, 2, ft.find(car_marker)), V)

bcu = [bcu_inflow, bcu_floor, bcu_walls, bcu_car]
bcp_outlet = dirichletbc(PETSc.ScalarType(0.0), locate_dofs_topological(Q, 2, ft.find(outlet_marker)), Q)
bcp = [bcp_outlet]

u = TrialFunction(V)
v = TestFunction(V)
u_ = Function(V)
u_.name = "Velocity"
u_s = Function(V)
u_n = Function(V)
p = TrialFunction(Q)
q = TestFunction(Q)
p_ = Function(Q)
p_.name = "Pressure"
p_n = Function(Q)
n = FacetNormal(mesh)
ds_tagged = ds(subdomain_data=ft)
f = Constant(mesh, np.zeros(3))

def epsilon(u):
    return sym(nabla_grad(u))

def sigma(u, p):
    return 2 * mu * epsilon(u) - p * Identity(3)


traction = dot(sigma(u_, p_), -n)
drag_form = form(traction[0] * ds_tagged(car_marker))
lift_form = form(traction[1] * ds_tagged(car_marker))

F1 = rho * dot((u - u_n) / k, v) * dx
F1 += rho * dot(dot(u_n, nabla_grad(u)), v) * dx
F1 += inner(sigma((u + u_n) / 2, p_n), epsilon(v)) * dx
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
L3 = form(dot(u_s, v) * dx - (k / rho) * dot(nabla_grad(p_ - p_n), v) * dx)
A3 = assemble_matrix(a3, bcs=bcu)
A3.assemble()
b3 = create_vector(extract_function_spaces(L3))


solver1 = PETSc.KSP().create(mesh.comm)
solver1.setOperators(A1)
solver1.setType("bcgs")
pc1 = solver1.getPC()
pc1.setType("hypre")
pc1.setHYPREType("boomeramg")
solver1.setTolerances(rtol=1e-5, atol=1e-8)

solver2 = PETSc.KSP().create(mesh.comm)
solver2.setOperators(A2)
solver2.setType("bcgs")
pc2 = solver2.getPC()
pc2.setType("hypre")
pc2.setHYPREType("boomeramg")
solver2.setTolerances(rtol=1e-8, atol=1e-12)

solver3 = PETSc.KSP().create(mesh.comm)
solver3.setOperators(A3)
solver3.setType("cg")
pc3 = solver3.getPC()
pc3.setType("sor")
solver3.setTolerances(rtol=1e-8, atol=1e-10)

# Расчёт
folder = Path(str(BASE_DIR / "results_benchmark"))
folder.mkdir(exist_ok=True, parents=True)
vtx_u = VTXWriter(mesh.comm, folder / "benchmark_u_2.bp", [u_], engine="BP4")
vtx_p = VTXWriter(mesh.comm, folder / "benchmark_p_2.bp", [p_], engine="BP4")
vtx_u.write(t)
vtx_p.write(t)

if mesh_comm.rank == 0:
    with open(folder / "benchmark_log_2.txt", "w") as f_log:
        f_log.write("Time, Drag, Lift, C_D, C_L, U_norm\n")

progress = tqdm.autonotebook.tqdm(desc="Solving PDE", total=num_steps)


for i in range(num_steps):
    if mesh_comm.rank == 0:
        progress.update(1)
    t += dt
    inlet_velocity.t = t
    u_inlet.interpolate(inlet_velocity)

    A1.zeroEntries()
    assemble_matrix(A1, a1, bcs=bcu)
    A1.assemble()
    solver1.setOperators(A1)

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

    drag_local = assemble_scalar(drag_form)
    drag_global = mesh.comm.allreduce(drag_local, op=MPI.SUM)

    lift_local = assemble_scalar(lift_form)
    lift_global = mesh.comm.allreduce(lift_local, op=MPI.SUM)


    C_D = (2.0 * drag_global) / (1.0 * U_mean**2 * A_ref)
    C_L = (2.0 * lift_global) / (1.0 * U_mean**2 * A_ref)

    if (i + 1) % save_freq == 0:
        u_norm = u_.x.petsc_vec.norm()

        vtx_u.write(t)
        vtx_p.write(t)

        if mesh.comm.rank == 0:
            with open(folder / "benchmark_log_2.txt", "a") as f_log:
                f_log.write(f"{t:.4f}, {drag_global:.6f}, {lift_global:.6f}, {C_D:.4f}, {C_L:.4f}, {u_norm:.6e}\n")

        r1 = solver1.getResidualNorm()
        r2 = solver2.getResidualNorm()
        it1 = solver1.getIterationNumber()
        it2 = solver2.getIterationNumber()

        if mesh.comm.rank == 0:
            print(f"t={t:.3f} | C_D={C_D:.4f} (ожидаем ~6.1) | C_L={C_L:.4f} (ожидаем ~0.009) | u_norm={u_norm:.3e}")

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