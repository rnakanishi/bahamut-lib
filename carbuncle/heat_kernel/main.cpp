#include <blas/grid_heat_kernel.h>
#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <iostream>

int main(int argc, char const *argv[]) {
  /* code */
  Ramuh::GridHeatKernel heat(
      Eigen::Array3i(25, 25, 1),
      Ramuh::BoundingBox3(Eigen::Array3d(0, 0, 0), Eigen::Array3d(1, 1, 0.1)));

  // initial conditions
  std::vector<Eigen::Array3i> position;
  std::vector<Eigen::Array3d> value;

  // position.emplace_back(4, 4, 0);
  // value.emplace_back(0, -1, 0);

  // position.emplace_back(5, 4, 0);
  // value.emplace_back(0, 1, 0);

  // position.emplace_back(4, 5, 0);
  // value.emplace_back(1, 0, 0);

  // position.emplace_back(4, 3, 0);
  // value.emplace_back(-1, 0, 0);

  // heat.addInitialCondition(heat.ijkToid(3, 3, 0), 0.2);

  // int id = 0;
  // for (auto p : position)
  //   heat.addInitialCondition(heat.ijkToid(p[0], p[1], p[2]), value[id++]);

  // auto result = heat.computeParallelTransport(10, false);

  // for (auto value : result)
  //   std::cerr << value.transpose() << std::endl;

  // heat.addInitialCondition(heat.ijkToid(12 + 1, 11 + 1, 0), 1);
  heat.addInitialCondition(heat.ijkToid(12 + 1, 12 + 1, 0), 1);
  // heat.addInitialCondition(heat.ijkToid(12 + 1, 13 + 1, 0), 1);
  // heat.addInitialCondition(heat.ijkToid(11, 12, 0), 1);
  // heat.addInitialCondition(heat.ijkToid(13, 12, 0), 1);

  // Find distance field
  auto distance = heat.computeHeatDistance();
  for (auto value : distance)
    std::cerr << value << " ";
  std::cerr << std::endl;

  return 0;
}