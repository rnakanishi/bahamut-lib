#include <structures/dual_cubes.h>
#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <utils/timer.hpp>
#include <sstream>
#include <iomanip>

int main(int argc, char const *argv[]) {
  Ramuh::DualCubes3 cubes(Eigen::Array3i(50, 50, 50),
                          Ramuh::BoundingBox3(Eigen::Array3d(-4., -2., -2.),
                                              Eigen::Array3d(4., 2., 2.)));
  cubes.initialize(Eigen::Array3d(-3, 0., 0), 0.5,
                   Ramuh::DualCubes3::ParametricSurface::CUBE);
  cubes.defineVelocity();
  // cubes.printCells();

  cubes.computeIntersection();
  cubes.computeNormals();
  cubes.extractSurface();
  Ramuh::Timer timer;
  for (int i = 1; i <= 100; i++) {
    timer.clearAll();
    timer.reset();
    cubes.cfl();
    do {
      cubes.advectWeno();
      // cubes.integrateLevelSet();
      timer.registerTime("Weno");
      cubes.computeIntersection();
      cubes.computeNormals();
      timer.registerTime("Hermite");
    } while (!cubes.advanceTime());
    cubes.extractSurface();
    timer.registerTime("Surface");
    timer.evaluateComponentsTime();
  }

  // auto surface = cubes.marchingTetrahedra();
  // Ramuh::FileWriter writer;
  // std::ostringstream objname;
  // objname << "results/marching/tetra_" << std::setfill('0') << std::setw(4)
  // << 0
  //         << ".obj";
  // writer.writeMeshModel(surface, objname.str());
  // }
  return 0;
}
