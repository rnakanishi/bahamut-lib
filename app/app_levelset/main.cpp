#include "levelset_fluid.h"
#include <iostream>
#include <structures/vector3.h>
#include <vector>
#include <cmath>

int main(void) {
  LevelSetFluid sim;
  sim.setResolution(Ramuh::Vector3i(64, 64, 64));

  std::cerr << "Initialized grid with size " << sim.domainSize()
            << ", resolution " << sim.resolution() << "  and h spacing "
            << sim.h() << std::endl;

  auto res = sim.resolution();
  auto h = sim.h();
  sim.addSphereSurface(Ramuh::Vector3d(0.3, 0.3, 0), 0.25);
  // for (int i = 0; i < res.x() + 1; i++)
  //   for (int j = 0; j < res.y() + 1; j++) {
  //     Ramuh::Vector3d pos(h * Ramuh::Vector3i(i, j, 0));
  //     sim[i][j][0] = std::cos(5.0 * pos.x()) * std::cos(5.0 * pos.y());
  //   }
  for (int i = 0; i < res.x() + 1; i++) {
    for (int j = 0; j < res.y() + 1; j++) {
      std::cerr << sim[i][j][0] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << "===\n";
  sim.setVelocity();
  sim.interpolateVelocitiesToVertices();
  // sim.printVertexVelocity();

  for (int i = 0; i < 5; i++) {
    sim.integrateLevelSet();
  }

  for (int i = 0; i < res.x() + 1; i++) {
    for (int j = 0; j < res.y() + 1; j++) {
      std::cerr << sim[i][j][0] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << "===\n";
  return 0;
}