#pragma once
#include <cmath>

namespace config {

int water_type = 0;
double water_m = 1.0;
double water_sigma = 1.0;
double water_epsilon = 1.0;

int salt_type = 1;
double salt_m = 5.0;
double salt_sigma = 2.0;
double salt_epsilon = 1.0;

double box_x = 20.0;
double box_y = 20.0;
double box_z = 20.0;
double mem_x = 2.0;

double T = 293.0;
double k_b = 1.0;
double water_mean_v = std::sqrt(k_b * T / water_m);
double salt_mean_v = std::sqrt(k_b * T / salt_m);
double concentration = 0.05;
int N = 130;
int salt_period = static_cast<int>(1 / concentration);
// int salt_N = static_cast<int>(N * 0.05);
// int water_N = N - salt_N;

} // namespace config
