#include <dolfin.h>
#include  <mshr.h>

#include "poisson.h"

using namespace dolfin;

// это источник слева снизу
class Source : public Expression {
  void eval(Array<double>& values, const Array<double>& x) const {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 100 / (dx * dx + dy * dy + 0.001);
  }
};

// это через верхнюю стенку условие на поток
class dUdN : public Expression {
  void eval(Array<double>& values, const Array<double>& x) const {
    values[0] = 1000 * (1 - abs(x[0] - 1.5));
  }
};

// это описание всех остальных стенок
// у них будет u = const
class DirichletBoundary : public SubDomain {
  bool inside(const Array<double>& x, bool on_boundary) const {
    return x[0] > 2 - DOLFIN_EPS || x[0] < DOLFIN_EPS || x[1] < DOLFIN_EPS || (
             abs(x[0] - 0.5) < 0.5 + DOLFIN_EPS && abs(x[1] - 1) < DOLFIN_EPS) || (
             abs(x[0] - 1) < DOLFIN_EPS && abs(x[1] - 1.5) < 0.5 + DOLFIN_EPS);
  }
};

int main() {
  auto big = mshr::Rectangle(Point(0.0, 0.0), Point(2.0, 2.0));
  auto small = mshr::Rectangle(Point(0.0, 1.0), Point(1.0, 2.0));
  auto dom = big - small;
  auto mesh = mshr::generate_mesh(dom, 256);
  auto V = std::make_shared<poisson::FunctionSpace>(mesh);

  auto u0 = std::make_shared<Constant>(-0.02);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  poisson::BilinearForm a(V, V);
  poisson::LinearForm L(V);
  auto f = std::make_shared<Source>();
  auto g = std::make_shared<dUdN>();
  L.f = f;
  L.g = g;

  Function u(V);
  solve(a == L, u, bc);

  File file("my_poisson.pvd");
  file << u;
}
