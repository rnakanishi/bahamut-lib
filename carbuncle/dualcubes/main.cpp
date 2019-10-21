#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <iomanip>
#include "dual_cubes.h"

int main(int argc, char const *argv[]) {
  Carbuncle::DualCubes3 cubes(
      Eigen::Array3i(40, 40, 40),
      Ramuh::BoundingBox3(Eigen::Array3d(-5, -5, -5), Eigen::Array3d(5, 5, 5)));
  cubes.initialize(Eigen::Array3d(0, 0, 0), 3,
                   Carbuncle::DualCubes3::ParametricSurface::CUBE);
  cubes.defineVelocity();
  // cubes.printCells();

  // cubes.computeIntersection();
  // cubes.analyticNormals(Eigen::Array3d(0, 0, 0), 0.5,
  //                       Carbuncle::DualCubes3::ParametricSurface::CUBE);
  // cubes.extractSurface();
  // return 1;

  // cubes.computeCellsGradient();
  // cubes.computeIntersectionAndNormals();
  // cubes.extractSurface();
  cubes.redistance();
  cubes.print();

  // return 1;

  for (int i = 1; i <= 15; i++) {
    cubes.advectWeno();
    // cubes.advectUpwind();

    if (i % 10 == 0)
      cubes.redistance();
    // cubes.computeIntersection();
    // cubes.computeNormals();

    cubes.computeWenoGradient();
    cubes.computeIntersectionAndNormals();
    cubes.extractSurface();
    cubes.print();
  }

  //   auto surface = cubes.marc hingTetrahedra();
  //   Ramuh::FileWriter writer;
  //   std::ostringstream objname;
  //   objname << "results/marching/tetra_" << std::setfill('0') << std::setw(4)
  //   << 0
  //   << ".obj";
  //   writer.writeMeshModel(surface, objname.str());
  // }
  return 0;
}
