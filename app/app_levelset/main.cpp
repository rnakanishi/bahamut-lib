#include "levelset_fluid.h"
#include <iostream>
#include <structures/vector3.h>
#include <vector>
#include <cmath>

int main(void) {
  LevelSetFluid sim;
  int resolution = 8;
  sim.setResolution(Ramuh::Vector3i(resolution, resolution, 1));
  sim.setSize(Ramuh::Vector3d(1.0, 1.0, 1.0 / resolution));

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

  sim.addGravity();
  sim.solvePressure();
  return 0;
}