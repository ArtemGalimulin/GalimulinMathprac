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
    if (type == config::water_type) {
      // Вода
      m = config::water_m;
      sigma = config::water_sigma;
      epsilon = config::water_epsilon;
    } else {
      // Соль
      m = config::salt_m;
      sigma = config::salt_sigma;
      epsilon = config::salt_epsilon;
    }
  }

  void move_1(double dt) {
    // Velocity Verlet part 1
    v += f * (dt / (2 * m));
    r += v * dt;
    f = {0.0, 0.0, 0.0};
  }

  void move_2(double dt) {
    // Velocity Verlet part 2
    v += f * (dt / (2 * m));
    f = {0.0, 0.0, 0.0};
  }
};

struct Box {
  double lx = config::box_x; // Полуширина по X
  double ly = config::box_y; // Полуширина по Y
  double lz = config::box_z; // Полуширина по Z

  bool inside(Vec3 r) const {
    return std::abs(r.x) < lx && std::abs(r.y) < ly && std::abs(r.z) < lz;
  }

  void reflect_particle(Particle& p) const {
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

  void reflect_particle(Particle& p) const {
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
      (type == config::water_type) ? config::water_mean_v : config::salt_mean_v;

  std::normal_distribution<double> dist(0.0, sigma);

  return {dist(gen), dist(gen), dist(gen)};
}

class Integrator {
  double current_mv2_left = 0.0;
  double current_mv2_right = 0.0;
  double current_virial_left = 0.0;
  double current_virial_right = 0.0;

  double accumulated_mv2_left = 0.0;
  double accumulated_mv2_right = 0.0;
  double accumulated_virial_left = 0.0;
  double accumulated_virial_right = 0.0;
  int steps_in_period = 0;

  int N_;
  int Nf_;
  double K_target_;
  double c_ = config::bussi_c;


  void calculate_forces(std::vector<Particle>& particles) {
    current_virial_left = 0.0;
    current_virial_right = 0.0;

    for (size_t i = 0; i < particles.size(); ++i) {
      for (size_t j = i + 1; j < particles.size(); ++j) {
        Particle& p1 = particles[i];
        Particle& p2 = particles[j];

        Vec3 r12 = p2.r - p1.r;
        double r2 = r12.norm_sq();
        if (r2 < config::singularity_threshold_2) {
          r2 = config::singularity_threshold_2;
        }

        // Параметры взаимодействия (смешивание по правилам Лоренца-Бертло)
        // sigma_ij = (sigma_i + sigma_j) / 2
        // epsilon_ij = sqrt(epsilon_i * epsilon_j)
        double sigma = (p1.sigma + p2.sigma) / 2.0;
        double epsilon = std::sqrt(p1.epsilon * p2.epsilon);

        double s2 = sigma * sigma;
        double r_inv2 = s2 / r2; // (sigma/r)^2
        double r_inv6 = r_inv2 * r_inv2 * r_inv2; // (sigma/r)^6

        // Сила Леннарда-Джонса: F = 24 * eps * [2*(sig/r)^12 - (sig/r)^6] / r^2
        // * r_vec
        double force_mod = 24.0 * epsilon * r_inv6 * (2.0 * r_inv6 - 1.0) / r2;
        Vec3 f12 = r12 * force_mod;
        p1.f -= f12;
        p2.f += f12;

        double v_pair = r12.dot(f12); // Вклад пары в вириал
        double mid_x = (p1.r.x + p2.r.x) / 2.0; // Средняя точка пары
        // Относим вириал к левому или правому отсеку
        if (mid_x < 0) {
          current_virial_left += v_pair;
        } else {
          current_virial_right += v_pair;
        }
      }
    }
  }

public:
  Integrator(int N, int Nf): N_(N), Nf_(Nf) {
    K_target_ = 0.5 * Nf_ * config::k_b * config::T;
  }

  void update(std::vector<Particle>& particles, Box& box, Membrane& membrane,
              double dt) {
    current_mv2_left = 0.0;
    current_mv2_right = 0.0;

    for (Particle& p : particles) {
      p.move_1(dt);
    }

    calculate_forces(particles);

    for (Particle& p : particles) {
      p.move_2(dt);
      box.reflect_particle(p);
      membrane.reflect_particle(p);

      if (p.r.x < 0) {
        current_mv2_left += p.m * p.v.norm_sq();
      } else {
        current_mv2_right += p.m * p.v.norm_sq();
      }
    }

    accumulated_mv2_left += current_mv2_left;
    accumulated_mv2_right += current_mv2_right;
    accumulated_virial_left += current_virial_left;
    accumulated_virial_right += current_virial_right;
    ++steps_in_period;
  }

  void apply_bussi_thermostat(std::vector<Particle>& particles) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<double> norm_dist(0.0, 1.0);
    static std::chi_squared_distribution<double> chi_dist(Nf_ - 1);

    // Формула А7 статьи Canonical sampling through velocity-rescaling
    double R1 = norm_dist(gen);
    double sum_Ri2 = chi_dist(gen);
    double K_old = 0.5 * (current_mv2_left + current_mv2_right);
    double B = K_target_ / (Nf_ * K_old);
    double A1 = B * (1 - c_) * (R1 * R1 + sum_Ri2);
    double A2 = 2 * std::sqrt(c_ * B * (1 - c_)) * R1;
    double alpha = std::sqrt(c_ + A1 + A2);

    for (auto& p : particles) {
      p.v *= alpha;
    }
  }

  void reset_accumulation() {
    accumulated_mv2_left = 0.0;
    accumulated_mv2_right = 0.0;
    accumulated_virial_left = 0.0;
    accumulated_virial_right = 0.0;
    steps_in_period = 0;
  }

  double get_T(int N) const {
    return (accumulated_mv2_left + accumulated_mv2_right) /
           (Nf_ * config::k_b * steps_in_period);
  }

  double get_left_P() const {
    if (steps_in_period == 0)
      return 0.0;
    // P = ( <2K> + <W> ) / 3V
    double P = (accumulated_mv2_left + accumulated_virial_left) /
               (3.0 * config::half_volume * steps_in_period);
    return P * config::P_to_bar; // Перевод в бары
  }

  double get_right_P() const {
    if (steps_in_period == 0)
      return 0.0;
    // P = ( <2K> + <W> ) / 3V
    double P = (accumulated_mv2_right + accumulated_virial_right) /
               (3.0 * config::half_volume * steps_in_period);
    return P * config::P_to_bar; // Перевод в бары
  }
};

class Simulation {
  std::string name_;
  Box box_;
  Membrane membrane_;
  std::vector<Particle> particles_;
  int N_;
  int Nf_;

  void rescale_velocities(double T_target) {
    double T_current = 0.0;
    for (const auto& p : particles_)
      T_current += p.m * p.v.norm_sq();
    T_current /= (Nf_ * config::k_b);

    double scale = std::sqrt(T_target / T_current);
    for (auto& p : particles_)
      p.v *= scale;
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

          int type = (count % config::salt_period == 0)
                       ? config::salt_type
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

      for (const auto& p : particles_) {
        out << p.type << " " << std::fixed << std::setprecision(6) << p.r.x
            << " " << p.r.y << " " << p.r.z << "\n";
      }

      out.close();
    } catch (const std::filesystem::filesystem_error& e) {
      std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
  }

public:
  Simulation(std::string& sim_name)
    : name_(sim_name), box_(Box()), membrane_(Membrane()), particles_({}) {
    particles_.reserve(config::N);
    setup_pure_water(config::N);
    setup_solution(config::N);
    N_ = static_cast<int>(particles_.size());
    Nf_ = 3 * N_ - 3;
    std::cout << "Эксперимент инициализирован. Итого частиц: " << N_
        << std::endl;
    rescale_velocities(config::T);
    std::cout << "Rescale выполнен, T -> " << config::T << " К\n";
  }

  void run() {
    Integrator integrator(N_, Nf_);

    snapshot(0);

    for (int step = 1; step <= config::total_steps; ++step) {
      integrator.update(particles_, box_, membrane_, config::dt);

      if (step % config::thermostat_period == 0) {
        integrator.apply_bussi_thermostat(particles_);
      }

      if (step % config::save_period == 0) {
        snapshot(step);
      }

      if (step % config::log_period == 0) {
        double p_l = integrator.get_left_P();
        double p_r = integrator.get_right_P();
        double t_mean = integrator.get_T(N_);

        std::cout << "Step: " << std::setw(6) << step << " / "
            << config::total_steps << " | T: " << std::fixed
            << std::setprecision(2) << t_mean << " K"
            << " | P_L: " << std::setw(8) << std::setprecision(2) << p_l
            << " bar"
            << " | P_R: " << std::setw(8) << std::setprecision(2) << p_r
            << " bar"
            << " | dP: " << std::setw(8) << (p_r - p_l) << " bar"
            << std::endl;

        // Сброс накопленных сумм для следующего периода усреднения
        integrator.reset_accumulation();
      }
    }
    std::cout << "Симуляция завершена успешно!" << std::endl;
  }
};

int main() {
  std::cout << "Введите название эксперимента:\n";
  std::string sim_name;
  std::cin >> sim_name;

  auto start_time = std::chrono::high_resolution_clock::now();
  Simulation simulation(sim_name);
  simulation.run();
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end_time - start_time;

  std::cout << "\n========================================\n";
  std::cout << "Общее время выполнения: " << std::fixed << std::setprecision(3)
      << duration.count() << " секунд." << std::endl;
  std::cout << "========================================\n";
}
