#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <iomanip>
#include "dual_cubes.h"

int main(int argc, char const *argv[]) {
  Carbuncle::DualCubes3 cubes(
      Eigen::Array3i(20, 20, 20),
      Ramuh::BoundingBox3(Eigen::Array3d(-2, -2, -2), Eigen::Array3d(2, 2, 2)));
  cubes.initialize(Eigen::Array3d(0, 0, 0), 0.5,
                   Carbuncle::DualCubes3::ParametricSurface::CUBE);
  cubes.defineVelocity();
  // cubes.printCells();

  cubes.computeIntersection();
  cubes.computeNormals();
  cubes.extractSurface();
  cubes.computeIntersection();
  cubes.print();
  return 1;
  for (int i = 1; i <= 1; i++) {
    cubes.advectSemiLagrangean();
    cubes.computeIntersection();
    cubes.computeNormals();
    cubes.extractSurface();
  }

  //   auto surface = cubes.marchingTetrahedra();
  //   Ramuh::FileWriter writer;
  //   std::ostringstream objname;
  //   objname << "results/marching/tetra_" << std::setfill('0') << std::setw(4)
  //   << 0
  //   << ".obj";
  //   writer.writeMeshModel(surface, objname.str());
  // }
  return 0;
}
