#pragma once
#include <cmath>

namespace config {
  // --- Временные параметры ---
  double dt = 0.001; // [пс] (пикосекунды) - шаг времени
  int total_steps = 10000; // [безразм] общее количество шагов
  int log_period = 2000; // [безразм] период вывода в консоль
  int save_period = 1000; // [безразм] период сохранения для OVITO
  int thermostat_period = 1000;

  double bussi_tau = 0.01;
  double bussi_c = std::exp(-dt * thermostat_period / bussi_tau);


  // --- Концентрация и количество ---
    int N = 500; // [безразм] общее число частиц
  int salt_period = 5; // [безразм] каждая такая частица - соль

  // --- Параметры ВОДЫ (H2O) ---
  int water_type = 0; // [безразм] идентификатор типа
  double water_m = 18.015; // [а.е.м.] (атомные единицы массы)
  double water_sigma = 3.166; // [Ангстрем] эффективный диаметр молекулы
  double water_epsilon = 0.650 * 100.0; // г*A^2/(моль*пс^2)

  // --- Параметры ЖЕЛТОЙ КРОВЯНОЙ СОЛИ (K4[Fe(CN)6]) ---
  int salt_type = 1; // [безразм] идентификатор типа
  double salt_m = 368.35; // [а.е.м.] масса всего комплекса
  double salt_sigma = 7.500; // [Ангстрем] диаметр аниона с оболочкой
  double salt_epsilon = 4.500 * 100.0; // г*A^2/(моль*пс^2)

  double singularity_threshold_2 =
      1.0; // [Ангстрем^2] порог для защиты от деления на 0

  // --- Геометрия системы ---
  double box_x = 80.0; // [Ангстрем] полудлина бокса по X
  double box_y = 40.0; // [Ангстрем] полуширина по Y
  double box_z = 40.0; // [Ангстрем] полувысота по Z
  double mem_x = 4.0; // [Ангстрем] полуширина мембраны
  double half_volume = 4.0 * box_x * box_y * box_z;
  double volume = 8.0 * box_x * box_y * box_z;

  // --- Термодинамика ---
  double T = 293.15; // [Кельвин] температура (20 градусов Цельсия)
  double k_b = 0.8314; // [г·Å²/(пс²·моль·К)] — R в единицах системы
  double water_mean_v =
      std::sqrt(k_b * T / water_m); // [Ангстрем/пс] средняя скорость
  double salt_mean_v =
      std::sqrt(k_b * T / salt_m); // [Ангстрем/пс] средняя скорость

  // Коэффициент перевода давления в бары для системы [г, A, пс]
  double P_to_bar = 166.054;
} // namespace config
