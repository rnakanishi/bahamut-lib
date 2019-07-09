#include <structures/dual_cubes.h>
#include <Eigen/Dense>

int main(int argc, char const *argv[]) {
  Ramuh::DualCubes3 cubes(16);
  cubes.initialize(Eigen::Array3d(0.5, 0.5, 0.5), 0.3);
  cubes.printCells();
  cubes.extractSurface();
  return 0;
}
