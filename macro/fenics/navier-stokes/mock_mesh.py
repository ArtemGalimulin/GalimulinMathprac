import gmsh
import numpy as np
import sys
from pathlib import Path

BASE_DIR = Path(__file__).parent


def generate_mock_msh(filename="pickup_mock.msh"):
    gmsh.initialize()
    gmsh.model.add("PickupMock")

    # Параметры канала
    L, W, H = 2.0, 1.0, 1.0
    # Параметры "машины" (коробочка внутри)
    cx0, cy0, cz0 = 0.5, 0.4, 0.0
    cl, cw, ch = 0.3, 0.2, 0.3

    # 1. Создаем внешний домен (воздух)
    air_box = gmsh.model.occ.addBox(0, 0, 0, L, W, H)

    # 2. Создаем "машину"
    car_box = gmsh.model.occ.addBox(cx0, cy0, cz0, cl, cw, ch)

    # 3. Вырезаем машину из воздушного домена
    # booleanFragments или cut для создания пустоты внутри
    gmsh.model.occ.cut([(3, air_box)], [(3, car_box)])
    gmsh.model.occ.synchronize()

    # 4. Определяем маркеры (Physical Groups)
    inlet_marker, outlet_marker, wall_marker, car_marker, air_marker = 1, 2, 3, 4, 5

    surfaces = gmsh.model.getEntities(2)

    inlet_list = []
    outlet_list = []
    wall_list = []
    car_list = []

    for surf in surfaces:
        tag = surf[1]
        # ВАЖНО: getCenterOfMass для OCC возвращает список [x, y, z]
        com = gmsh.model.occ.getCenterOfMass(2, tag)

        x, y, z = com[0], com[1], com[2]

        if np.isclose(x, 0):
            inlet_list.append(tag)
        elif np.isclose(x, L):
            outlet_list.append(tag)
        elif np.isclose(y, 0) or np.isclose(y, W) or np.isclose(z, H) or np.isclose(z, 0):
            # Проверяем, не является ли эта поверхность частью "машины"
            # Если Z=0, но координаты X,Y внутри границ машины — это дно машины
            if np.isclose(z, 0) and (cx0 <= x <= cx0 + cl) and (cy0 <= y <= cy0 + cw):
                car_list.append(tag)
            else:
                wall_list.append(tag)
        else:
            car_list.append(tag)

    gmsh.model.addPhysicalGroup(2, inlet_list, inlet_marker, name="Inlet")
    gmsh.model.addPhysicalGroup(2, outlet_list, outlet_marker, name="Outlet")
    gmsh.model.addPhysicalGroup(2, wall_list, wall_marker, name="Walls")
    gmsh.model.addPhysicalGroup(2, car_list, car_marker, name="Car")

    # Объем
    volumes = gmsh.model.getEntities(3)
    gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes], air_marker, name="Air")

    # 5. Сетка
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)

    gmsh.model.mesh.generate(3)
    gmsh.write(str(BASE_DIR / filename))

    print(f"Файл {filename} успешно создан!")
    gmsh.finalize()


if __name__ == "__main__":
    generate_mock_msh()
