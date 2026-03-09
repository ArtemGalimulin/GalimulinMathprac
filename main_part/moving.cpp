#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

class CalcNode {
  friend class CalcMesh;

protected:
  double x;
  double y;
  double z;
  double x0;
  double y0;
  double z0;
  double smth;
  double vx;
  double vy;
  double vz;

public:
  CalcNode() : x(0.0), y(0.0), z(0.0),
               x0(0.0), y0(0.0), z0(0.0),
               smth(0.0), vx(0.0), vy(0.0), vz(0.0) {
  }

  CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz)
    : x(x), y(y), z(z), x0(x), y0(y), z0(z), smth(smth), vx(vx), vy(vy), vz(vz) {
  }

  void move(double tau, unsigned int step) {
    // белка машет хвостом вдали от точки его присоединения как
    // вокруг оси, а вблизи перехода там плавный переход, задающийся весом
    // поэтому хвост не как отдельное, а будто присоединен к белке
    // этот вес как раз и отображается цветом
    double t = step * tau;
    double omega = 3.0;
    double A = 0.5;
    double dist_to_plane = (z0 - (53.0 - 3.6 * y0)) / sqrt(1.0 + 3.6 * 3.6);
    double transition = 6.0;
    double weight = 0.5 * (1.0 + tanh(dist_to_plane / transition));
    double theta = weight * A * sin(omega * t);
    double dtheta = weight * A * omega * cos(omega * t);
    double dist = sqrt(x0 * x0 + (y0 - 15.0) * (y0 - 15.0));
    double phi0 = atan2(y0 - 15.0, x0);

    x = dist * cos(phi0 + theta);
    y = dist * sin(phi0 + theta) + 15.0;
    vx = dist * dtheta * (-sin(phi0 + theta));
    vy = dist * dtheta * (cos(phi0 + theta));

    smth = weight;
  }
};

class Element {
  friend class CalcMesh;

protected:
  unsigned long nodesIds[4];
};


class CalcMesh {
protected:
  vector<CalcNode> nodes;
  vector<Element> elements;

public:
  CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {
    nodes.resize(nodesCoords.size() / 3);
    for (unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
      double pointX = nodesCoords[i * 3];
      double pointY = nodesCoords[i * 3 + 1];
      double pointZ = nodesCoords[i * 3 + 2];
      // double smth = sin(pointX * 0.1);
      double smth = 1;
      nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
    }

    elements.resize(tetrsPoints.size() / 4);
    for (unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
      elements[i].nodesIds[0] = tetrsPoints[i * 4] - 1;
      elements[i].nodesIds[1] = tetrsPoints[i * 4 + 1] - 1;
      elements[i].nodesIds[2] = tetrsPoints[i * 4 + 2] - 1;
      elements[i].nodesIds[3] = tetrsPoints[i * 4 + 3] - 1;
    }
  }

  void doTimeStep(double tau, unsigned int step) {
    for (unsigned int i = 0; i < nodes.size(); ++i) {
      nodes[i].move(tau, step);
    }
  }

  void snapshot(unsigned int snap_number) {
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    auto smth = vtkSmartPointer<vtkDoubleArray>::New();
    smth->SetName("tail measure");

    auto vel = vtkSmartPointer<vtkDoubleArray>::New();
    vel->SetName("velocity");
    vel->SetNumberOfComponents(3);

    for (unsigned int i = 0; i < nodes.size(); i++) {
      dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

      double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
      vel->InsertNextTuple(_vel);

      smth->InsertNextValue(nodes[i].smth);
    }

    unstructuredGrid->SetPoints(dumpPoints);
    unstructuredGrid->GetPointData()->AddArray(vel);
    unstructuredGrid->GetPointData()->AddArray(smth);

    for (unsigned int i = 0; i < elements.size(); i++) {
      auto tetra = vtkSmartPointer<vtkTetra>::New();
      tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
      tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
      tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
      tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
      unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
    }

    string fileName = "moving-" + std::to_string(snap_number) + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetDataModeToAscii(); // это важно, а то иначе 3д не работает
    writer->SetInputData(unstructuredGrid);
    writer->Write();
  }
};

int main() {
  double h = 6.0;
  double tau = 0.01;

  const unsigned int GMSH_TETR_CODE = 4;

  gmsh::initialize();
  gmsh::model::add("squirrel");

  try {
    gmsh::merge("../squirrel.stl");
  } catch (...) {
    gmsh::logger::write("Could not load STL mesh: bye!");
    gmsh::finalize();
    return -1;
  }

  double angle = 40;
  bool forceParametrizablePatches = true;
  bool includeBoundary = true;
  double curveAngle = 180;
  gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                      forceParametrizablePatches,
                                      curveAngle * M_PI / 180.);
  gmsh::model::mesh::createGeometry();

  std::vector<std::pair<int, int>> s;
  gmsh::model::getEntities(s, 2);
  std::vector<int> sl;
  for (auto surf: s) sl.push_back(surf.second);
  int l = gmsh::model::geo::addSurfaceLoop(sl);
  gmsh::model::geo::addVolume({l});

  gmsh::model::geo::synchronize();

  int f = gmsh::model::mesh::field::add("MathEval");
  gmsh::model::mesh::field::setString(f, "F", std::to_string(h));
  gmsh::model::mesh::field::setAsBackgroundMesh(f);

  gmsh::model::mesh::generate(3);

  std::vector<double> nodesCoord;
  std::vector<std::size_t> nodeTags;
  std::vector<double> parametricCoord;
  gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

  std::vector<std::size_t>* tetrsNodesTags = nullptr;
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t>> elementTags;
  std::vector<std::vector<std::size_t>> elementNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
  for (unsigned int i = 0; i < elementTypes.size(); i++) {
    if (elementTypes[i] != GMSH_TETR_CODE)
      continue;
    tetrsNodesTags = &elementNodeTags[i];
  }

  if (tetrsNodesTags == nullptr) {
    cout << "Can not find tetra data. Exiting." << endl;
    gmsh::finalize();
    return -2;
  }

  cout << "The model has " << nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

  for (int i = 0; i < nodeTags.size(); ++i) {
    assert(i == nodeTags[i] - 1);
  }
  assert(tetrsNodesTags->size() % 4 == 0);

  CalcMesh mesh(nodesCoord, *tetrsNodesTags);

  // gmsh::finalize();

  for (unsigned int step = 1; step < 220; ++step) {
    mesh.doTimeStep(tau, step);
    mesh.snapshot(step);
  }

  return 0;
}
