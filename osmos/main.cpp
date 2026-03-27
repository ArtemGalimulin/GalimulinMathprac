#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "Particle.cpp"
#include "Vec3.cpp"
#include "config.cpp"

struct Particle {
  Vec3 r, v, f;
  double m, sigma, epsilon;
  int type;

  Particle(double x, double y, double z, int t) : r(x, y, z), type(t) {
    if (type == config::water_type) { // Вода
      m = config::water_m;
      sigma = config::water_sigma;
      epsilon = config::water_epsilon;
    } else { // Соль
      m = config::salt_m;
      sigma = config::salt_sigma;
      epsilon = config::salt_epsilon;
    }
  }

  void move(double dt) {
    v += f * (dt / m);
    r += v * dt;
    f = Vec3(0, 0, 0);
  }
};

struct Box {
  double lx = config::box_x; // Полуширина по X
  double ly = config::box_y; // Полуширина по Y
  double lz = config::box_z; // Полуширина по Z

  bool inside(Vec3 r) const {
    return std::abs(r.x) < lx && std::abs(r.y) < ly && std::abs(r.z) < lz;
  }

  void reflect_particle(Particle &p) const {
    if (std::abs(p.r.x) >= lx) {
      p.v.x *= -1.0;
      p.r.x = (p.r.x > 0) ? (2.0 * lx - p.r.x) : (-2.0 * lx - p.r.x);
    }

    if (std::abs(p.r.y) >= ly) {
      p.v.y *= -1.0;
      p.r.y = (p.r.y > 0) ? (2.0 * ly - p.r.y) : (-2.0 * ly - p.r.y);
    }

    if (std::abs(p.r.z) >= lz) {
      p.v.z *= -1.0;
      p.r.z = (p.r.z > 0) ? (2.0 * lz - p.r.z) : (-2.0 * lz - p.r.z);
    }
  }
};

struct Membrane {
  double hw = config::mem_x;

  bool is_inside(Vec3 r) const { return std::abs(r.x) < hw; }

  void reflect_particle(Particle &p) const {
    if (p.type == config::salt_type) {
      if (is_inside(p.r)) {
        p.v.x *= -1.0;
        p.r.x = (p.r.x > 0) ? (2.0 * hw - p.r.x) : (-2.0 * hw - p.r.x);
      }
    }
  }
};

Vec3 generate_maxwell_velocity(int type) {
  static std::random_device rd;
  static std::mt19937 gen(rd());

  double sigma =
      (type == config::water_type) ? config::water_sigma : config::salt_sigma;

  std::normal_distribution<double> dist(0.0, sigma);

  return {dist(gen), dist(gen), dist(gen)};
}

class Simulation {
  Box box_;
  Membrane membrane_;
  std::vector<Particle> particles_;

  void setup(int N, int concentration) {
    double lx = box_.lx - membrane_.hw;
    double ly = box_.ly;
    double lz = box_.lz;
    double coef =
        std::pow(static_cast<double>(N) / 2 / (lx * ly * lz), 1.0 / 3.0);
    int nx = static_cast<int>(coef * lx);
  }

  void setup_pure_water(int N) {
    double lx_eff = box_.lx - membrane_.hw;
    double ly_eff = 2.0 * box_.ly;
    double lz_eff = 2.0 * box_.lz;

    double coef = std::pow(
        static_cast<double>(N) / 2.0 / (lx_eff * ly_eff * lz_eff), 1.0 / 3.0);

    int nx = static_cast<int>(coef * lx_eff);
    int ny = static_cast<int>(coef * ly_eff);
    int nz = static_cast<int>(coef * lz_eff);

    double dx = lx_eff / nx;
    double dy = ly_eff / ny;
    double dz = lz_eff / nz;

    for (int i = 0; i < nx; ++i) {
      double x = -box_.lx + (i + 0.5) * dx;
      for (int j = 0; j < ny; ++j) {
        double y = -box_.ly + (j + 0.5) * dy;
        for (int k = 0; k < nz; ++k) {
          double z = -box_.lz + (k + 0.5) * dz;

          particles_.emplace_back(x, y, z, config::water_type);
          particles_.back().v = generate_maxwell_velocity(config::water_type);
        }
      }
    }
  }

  void setup_solution(int N) {
    double lx_eff = box_.lx - membrane_.hw;
    double ly_eff = 2.0 * box_.ly;
    double lz_eff = 2.0 * box_.lz;

    double coef = std::pow(
        static_cast<double>(N) / 2.0 / (lx_eff * ly_eff * lz_eff), 1.0 / 3.0);

    int nx = static_cast<int>(coef * lx_eff);
    int ny = static_cast<int>(coef * ly_eff);
    int nz = static_cast<int>(coef * lz_eff);

    double dx = lx_eff / nx;
    double dy = ly_eff / ny;
    double dz = lz_eff / nz;

    int count = 0;
    for (int i = 0; i < nx; ++i) {
      double x = membrane_.hw + (i + 0.5) * dx;
      for (int j = 0; j < ny; ++j) {
        double y = -box_.ly + (j + 0.5) * dy;
        for (int k = 0; k < nz; ++k) {
          double z = -box_.lz + (k + 0.5) * dz;

          int type = (count % config::salt_period == 0) ? config::salt_type
                                                        : config::water_type;

          particles_.emplace_back(x, y, z, type);
          particles_.back().v = generate_maxwell_velocity(type);
          ++count;
        }
      }
    }
  }

public:
  void save_snapshot(int step) {
    // Меняем расширение на .xyz
    std::string filename = "snap_" + std::to_string(step) + ".xyz";
    std::ofstream out(filename);

    if (!out.is_open()) {
      std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
      return;
    }

    // 1. Первая строка: количество частиц
    out << particles_.size() << "\n";

    // 2. Вторая строка: комментарий (номер шага симуляции)
    out << "Lattice snapshot at step " << step << "\n";

    // 3. Остальные строки: Type X Y Z (через пробел)
    for (const Particle &p : particles_) {
      // В формате XYZ тип частицы обычно идет первым (как символ элемента)
      // Но OVITO поймет и числовой тип.
      out << p.type << " "
          << std::fixed << std::setprecision(6)
          << p.r.x << " " << p.r.y << " " << p.r.z << "\n";
    }

    out.close();
  }

public:
  Simulation() : box_(Box()), membrane_(Membrane()), particles_({}) {
    particles_.reserve(config::N);

    std::cout << "Симуляция создана. Размеры коробки: " << box_.lx * 2 << "x"
              << box_.ly * 2 << "x" << box_.lz * 2 << std::endl;
  }

  void init_experiment(int target_N) {
    particles_.clear();
    setup_pure_water(target_N);
    setup_solution(target_N);
    std::cout << "Эксперимент инициализирован. Итого частиц: "
              << particles_.size() << std::endl;
  }

  const std::vector<Particle> &get_particles() const { return particles_; }
};

int main() {
  // std::cout << 123;

  Simulation simulation;
  simulation.init_experiment(config::N);
  simulation.save_snapshot(0);
}