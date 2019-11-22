#include <blas/grid_heat_kernel.h>
#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <iostream>

int main(int argc, char const *argv[]) {
  /* code */
  Ramuh::GridHeatKernel heat(Eigen::Array3i(50), Ramuh::BoundingBox3(0, 1));

  // initial conditions
  std::vector<Eigen::Array3i> position;
  std::vector<Eigen::Array3d> value;

  position.emplace_back(46, 49, 49);
  value.emplace_back(1, 1, 1);

  position.emplace_back(45, 49, 49);
  value.emplace_back(-1, -1, -1);

  int id = 0;
  for (auto p : position)
    heat.addInitialCondition(heat.ijkToid(p[0], p[1], p[2]), value[id++]);

  auto result = heat.solve();

  for (auto value : result)
    std::cerr << value.transpose() << std::endl;

  return 0;
}
