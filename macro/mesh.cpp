#include <gmsh.h>
#include <vector>
#include <iostream>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("DetailedPickup");

    // --- ПАРАМЕТРЫ (в метрах) ---
    double W = 2.0;               // Общая ширина
    double L_hood = 1.2;          // Длина капота
    double L_cabin = 2.0;         // Длина кабины
    double L_bed = 2.2;           // Длина кузова

    double H_ground = 0.4;        // Клиренс
    double H_hood = 1.0;          // Высота капота (от рамы)
    double H_cabin = 1.8;         // Высота кабины (от рамы)
    double H_bed = 1.0;           // Высота бортов (от рамы)

    double R_wheel = 0.45;        // Радиус колес
    double W_wheel = 0.3;         // Ширина колес

    namespace occ = gmsh::model::occ;

    // 1. СОЗДАНИЕ ОСНОВНЫХ ОБЪЕМОВ (с небольшим перекрытием 0.01 для стабильности Boolean)
    // Капот
    int hood = occ::addBox(0, 0, H_ground, L_hood + 0.01, W, H_hood);

    // Кабина
    int cabin = occ::addBox(L_hood, 0, H_ground, L_cabin + 0.01, W, H_cabin);

    // Багажник (внешний)
    int bed = occ::addBox(L_hood + L_cabin, 0, H_ground, L_bed, W, H_bed);

    // 2. ОБЪЕДИНЕНИЕ КУЗОВА
    std::vector<std::pair<int, int>> body;
    std::vector<std::vector<std::pair<int, int>>> map;
    occ::fuse({{3, hood}}, {{3, cabin}, {3, bed}}, body, map);

    // 3. СОЗДАНИЕ ВЫРЕЗОВ (CUTTERS)
    std::vector<std::pair<int, int>> cutters;

    // Скос лобового стекла (наклонный бокс)
    int glass_cutter = occ::addBox(L_hood - 0.2, -0.1, H_ground + H_hood + 0.1, 1.5, W + 0.2, 1.0);
    occ::rotate({{3, glass_cutter}}, L_hood, 0, H_ground + H_hood, 0, 1, 0, -0.5); // Наклон 0.5 рад
    cutters.push_back({3, glass_cutter});

    // Полость багажника
    double wall = 0.1;
    int bed_inner = occ::addBox(L_hood + L_cabin + wall, wall, H_ground + wall, L_bed - wall*2, W - wall*2, H_bed + 0.5);
    cutters.push_back({3, bed_inner});

    // Колесные арки (цилиндры по бокам)
    double wheel_track = 0.2; // На сколько колеса утоплены
    double wheel_base_front = L_hood * 0.6;
    double wheel_base_rear = L_hood + L_cabin + L_bed * 0.7;

    auto add_arch = [&](double x, double y_side) {
        int cyl = occ::addCylinder(x, y_side, H_ground, 0, (y_side > 1 ? -0.5 : 0.5), 0, R_wheel + 0.1);
        cutters.push_back({3, cyl});
    };

    add_arch(wheel_base_front, -0.1);    // Перед лево
    add_arch(wheel_base_front, W + 0.1); // Перед право
    add_arch(wheel_base_rear, -0.1);     // Зад лево
    add_arch(wheel_base_rear, W + 0.1);  // Зад право

    // 4. ПРИМЕНЕНИЕ ВЫРЕЗОВ
    std::vector<std::pair<int, int>> final_body;
    occ::cut(body, cutters, final_body, map);

    // 5. КОЛЕСА (отдельные объекты, не впаянные в кузов для чистоты сетки)
    auto add_wheel = [&](double x, double y) {
        occ::addCylinder(x, y, H_ground, 0, (y > 1 ? -W_wheel : W_wheel), 0, R_wheel);
    };
    add_wheel(wheel_base_front, 0);
    add_wheel(wheel_base_front, W);
    add_wheel(wheel_base_rear, 0);
    add_wheel(wheel_base_rear, W);

    occ::synchronize();

    // 6. СЕТКА (настройка размеров)
    double lc = 0.12;
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", lc * 0.5);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", lc);

    gmsh::model::mesh::generate(3);
    gmsh::write("pickup_v2.msh");

    // Запуск GUI
    bool nopopup = false;
    for(int i=0; i<argc; i++) if(std::string(argv[i]) == "-nopopup") nopopup = true;
    if(!nopopup) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}