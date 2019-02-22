#include "levelset_fluid.h"
#include <iostream>
#include <structures/vector3.h>
#include <vector>
#include <cmath>

int main(void) {
  LevelSetFluid sim;
  sim.setResolution(Ramuh::Vector3i(4, 4, 1));

  std::cerr << "Initialized grid with size " << sim.resolution()
            << "  and resolution " << sim.h() << std::endl;

  auto res = sim.resolution();
  auto h = sim.h();
  for (int i = 0; i < res.x() + 1; i++)
    for (int j = 0; j < res.y() + 1; j++) {
      Ramuh::Vector3d pos(h * Ramuh::Vector3i(i, j, 0));
      sim[i][j][0] = std::cos(5.0 * pos.x()) * std::cos(5.0 * pos.y());
    }
  for (int i = 0; i < res.x() + 1; i++) {
    for (int j = 0; j < res.y() + 1; j++) {
      std::cerr << sim[i][j][0] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << "===\n";
  sim.integrateLevelSet();

  sim.setVelocity();
  sim.printFaceVelocity();
  sim.interpolateVelocitiesToVertices();
  sim.printVertexVelocity();

  return 0;
}