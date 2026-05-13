#include <set>
#include <cmath>
#include <gmsh.h>

int main(int argc, char** argv) {
  gmsh::initialize();
  gmsh::model::add("Model");

  try {
    gmsh::merge("../dodge_light.stl");
  } catch (...) {
    gmsh::logger::write("Could not load STL mesh: bye!");
    gmsh::finalize();
    return 0;
  }

  gmsh::model::mesh::classifySurfaces(80 * M_PI / 180., true, true,
                                      180 * M_PI / 180.);
  gmsh::model::mesh::createGeometry();

  std::vector<std::pair<int, int>> s;
  gmsh::model::getEntities(s, 2);
  std::vector<int> sl;
  for (auto surf : s) sl.push_back(surf.second);
  int carLoop = gmsh::model::geo::addSurfaceLoop(sl);
  int volCar = gmsh::model::geo::addVolume({carLoop});

  // Машина: X[-2.32, 0.43]  Y[-2.56, 3.25]  Z[0.0, 2.17]
  double xmin = -2.32 - 2;
  double xmax = 0.43 + 2;
  double ymin = -2.56 - 3.0;
  double ymax = 3.25 + 5.0;
  double zmin = -0.5;
  double zmax = 2.17 + 2;

  int p1 = gmsh::model::geo::addPoint(xmin, ymin, zmin);
  int p2 = gmsh::model::geo::addPoint(xmax, ymin, zmin);
  int p3 = gmsh::model::geo::addPoint(xmax, ymax, zmin);
  int p4 = gmsh::model::geo::addPoint(xmin, ymax, zmin);
  int p5 = gmsh::model::geo::addPoint(xmin, ymin, zmax);
  int p6 = gmsh::model::geo::addPoint(xmax, ymin, zmax);
  int p7 = gmsh::model::geo::addPoint(xmax, ymax, zmax);
  int p8 = gmsh::model::geo::addPoint(xmin, ymax, zmax);

  int l1 = gmsh::model::geo::addLine(p1, p2);
  int l2 = gmsh::model::geo::addLine(p2, p3);
  int l3 = gmsh::model::geo::addLine(p3, p4);
  int l4 = gmsh::model::geo::addLine(p4, p1);
  int l5 = gmsh::model::geo::addLine(p5, p6);
  int l6 = gmsh::model::geo::addLine(p6, p7);
  int l7 = gmsh::model::geo::addLine(p7, p8);
  int l8 = gmsh::model::geo::addLine(p8, p5);
  int l9 = gmsh::model::geo::addLine(p1, p5);
  int l10 = gmsh::model::geo::addLine(p2, p6);
  int l11 = gmsh::model::geo::addLine(p3, p7);
  int l12 = gmsh::model::geo::addLine(p4, p8);

  int cl1 = gmsh::model::geo::addCurveLoop({l1, l2, l3, l4});
  int cl2 = gmsh::model::geo::addCurveLoop({l5, l6, l7, l8});
  int cl3 = gmsh::model::geo::addCurveLoop({l1, l10, -l5, -l9});
  int cl4 = gmsh::model::geo::addCurveLoop({l2, l11, -l6, -l10});
  int cl5 = gmsh::model::geo::addCurveLoop({l3, l12, -l7, -l11});
  int cl6 = gmsh::model::geo::addCurveLoop({l4, l9, -l8, -l12});

  int s1 = gmsh::model::geo::addPlaneSurface({cl1});
  int s2 = gmsh::model::geo::addPlaneSurface({cl2});
  int s3 = gmsh::model::geo::addPlaneSurface({cl3});
  int s4 = gmsh::model::geo::addPlaneSurface({cl4});
  int s5 = gmsh::model::geo::addPlaneSurface({cl5});
  int s6 = gmsh::model::geo::addPlaneSurface({cl6});

  int tunnelLoop = gmsh::model::geo::addSurfaceLoop({s1, s2, s3, s4, s5, s6});
  int volAir = gmsh::model::geo::addVolume({tunnelLoop, carLoop});

  gmsh::model::geo::synchronize();

  gmsh::model::addPhysicalGroup(2, {s3}, 1, "inlet");
  gmsh::model::addPhysicalGroup(2, {s5}, 2, "outlet");
  gmsh::model::addPhysicalGroup(2, {s2, s4, s6}, 3, "walls");
  gmsh::model::addPhysicalGroup(2, {s1}, 4, "floor");
  gmsh::model::addPhysicalGroup(2, sl, 5, "car_wall");
  gmsh::model::addPhysicalGroup(3, {volAir}, 6, "air");

  // int f = gmsh::model::mesh::field::add("MathEval");
  // gmsh::model::mesh::field::setString(f, "F", "0.3");
  // gmsh::model::mesh::field::setAsBackgroundMesh(f);

  int f_dist = gmsh::model::mesh::field::add("Distance");
  gmsh::model::mesh::field::setNumbers(f_dist, "SurfacesList", std::vector<double>(sl.begin(), sl.end()));

  int f_thresh = gmsh::model::mesh::field::add("Threshold");
  gmsh::model::mesh::field::setNumber(f_thresh, "InField", f_dist);
  gmsh::model::mesh::field::setNumber(f_thresh, "SizeMin", 0.08); // Оптимально для скорости
  gmsh::model::mesh::field::setNumber(f_thresh, "SizeMax", 1.0);
  gmsh::model::mesh::field::setNumber(f_thresh, "DistMin", 0.2);
  gmsh::model::mesh::field::setNumber(f_thresh, "DistMax", 1.5);

  gmsh::model::mesh::field::setAsBackgroundMesh(f_thresh);

  // Эти опции обязательны, чтобы Gmsh не игнорировал твои настройки
  gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
  gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

  gmsh::option::setNumber("Mesh.MaxNumThreads3D", 8);
  gmsh::option::setNumber("Mesh.Algorithm3D", 1);
  gmsh::option::setNumber("Mesh.Optimize", 0);

  gmsh::model::mesh::generate(3);
  gmsh::write("light_03_relative.msh");

  std::set<std::string> args(argv, argv + argc);
  if (!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();
  return 0;
}
