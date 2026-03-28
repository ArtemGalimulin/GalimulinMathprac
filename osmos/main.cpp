#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

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

class Integrator {
  void calculate_forces(std::vector<Particle> &particles) {
    for (Particle &p : particles)
      p.f = Vec3(0, 0, 0);

    for (size_t i = 0; i < particles.size(); ++i) {
      for (size_t j = i + 1; j < particles.size(); ++j) {
        Particle &p1 = particles[i];
        Particle &p2 = particles[j];

        Vec3 r12 = p2.r - p1.r;
        double r2 = r12.norm_sq();
        if (r2 < config::singularity_threshold) {
          r2 = config::singularity_threshold;
        }

        // Параметры взаимодействия (смешивание по правилам Лоренца-Бертло)
        // sigma_ij = (sigma_i + sigma_j) / 2
        // epsilon_ij = sqrt(epsilon_i * epsilon_j)
        double sigma = (p1.sigma + p2.sigma) / 2.0;
        double epsilon = std::sqrt(p1.epsilon * p2.epsilon);

        double s2 = sigma * sigma;
        double r_inv2 = s2 / r2;                  // (sigma/r)^2
        double r_inv6 = r_inv2 * r_inv2 * r_inv2; // (sigma/r)^6

        // Сила Леннарда-Джонса: F = 24 * eps * [2*(sig/r)^12 - (sig/r)^6] / r^2
        // * r_vec
        double force_mod = 24.0 * epsilon * r_inv6 * (2.0 * r_inv6 - 1.0) / r2;
        Vec3 f12 = r12 * force_mod;
        p1.f -= f12;
        p2.f += f12;
      }
    }
  }

public:
  void update(std::vector<Particle> &particles, Box &box, Membrane &membrane,
              double dt) {
    calculate_forces(particles);
    for (Particle &p : particles) {
      p.move(dt);
      box.reflect_particle(p);
      membrane.reflect_particle(p);
    }
  }
};

class Simulation {
  std::string name_;
  Box box_;
  Membrane membrane_;
  std::vector<Particle> particles_;

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

  void snapshot(int step) {
    try {
      std::filesystem::path dir(name_);
      if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
      }

      std::filesystem::path filepath =
          dir / ("snap_" + std::to_string(step) + ".xyz");

      std::ofstream out(filepath);
      if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filepath
                  << std::endl;
        return;
      }

      out << particles_.size() << "\n";
      out << "Simulation: " << name_ << " | Step: " << step << "\n";

      for (const auto &p : particles_) {
        out << p.type << " " << std::fixed << std::setprecision(6) << p.r.x
            << " " << p.r.y << " " << p.r.z << "\n";
      }

      out.close();

    } catch (const std::filesystem::filesystem_error &e) {
      std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
  }

public:
  Simulation(std::string &sim_name)
      : name_(sim_name), box_(Box()), membrane_(Membrane()), particles_({}) {
    particles_.reserve(config::N);
    setup_pure_water(config::N);
    setup_solution(config::N);
    std::cout << "Эксперимент инициализирован. Итого частиц: "
              << particles_.size() << std::endl;
  }

  void run() {
    Integrator integrator;

    snapshot(0);

    for (int step = 1; step <= config::total_steps; ++step) {
      integrator.update(particles_, box_, membrane_, config::dt);

      if (step % config::log_period == 0) {
        std::cout << "Step: " << step << " / " << config::total_steps
                  << std::endl;
      }

      if (step % config::save_period == 0) {
        snapshot(step);
      }
    }
    std::cout << "Симуляция завершена успешно!" << std::endl;
  }
};

int main() {
  std::cout << "Введите название эксперимента:\n";
  std::string sim_name;
  std::cin >> sim_name;

  Simulation simulation(sim_name);
  simulation.run();
}