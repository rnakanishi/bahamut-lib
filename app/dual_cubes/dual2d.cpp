#include <structures/dual_cubes.h>
#include <geometry/bounding_box.h>
#include <Eigen/Dense>

int main(int argc, char const *argv[]) {
  Ramuh::DualCubes3 cubes(
      Eigen::Array3i(16, 16, 16),
      Ramuh::BoundingBox3(Eigen::Array3d(-1, -1, -1), Eigen::Array3d(1, 1, 1)));
  cubes.initialize(Eigen::Array3d(0, 0, 0), 0.3);
  cubes.printCells();
  cubes.extractSurface();
  return 0;
}
