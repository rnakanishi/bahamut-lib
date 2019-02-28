#include "levelset_fluid.h"
#include <iostream>
#include <geometry/vector3.h>
#include <geometry/matrix.h>
#include <vector>
#include <cmath>
#include <utils/file_writer.h>

int main(void) {
  LevelSetFluid sim;
  Ramuh::FileWriter output;
  int resolution = 8;
  sim.setResolution(Ramuh::Vector3i(resolution, resolution, 1));
  sim.setSize(Ramuh::Vector3d(1.0, 1.0, 1.0 / resolution));

  std::cerr << "Initialized grid with size " << sim.domainSize()
            << ", resolution " << sim.resolution() << "  and h spacing "
            << sim.h() << std::endl;

  auto res = sim.resolution();
  auto h = sim.h();
  sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.1, 0), 0.3);

  for (int j = 0; j < res.y(); j++) {
    for (int i = 0; i < res.x(); i++) {
      std::cerr << sim[i][j][0] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << "===\n";
  sim.setVelocity();
  // sim.printVertexVelocity();

  for (int frame = 0; frame < 3; frame++) {
    sim.advectGridVelocity();
    sim.addGravity();
    sim.boundaryVelocities();
    sim.solvePressure();
    sim.interpolateVelocitiesToVertices();
    sim.integrateLevelSet();

    output.writeLevelSet(sim, "data/" + frame);
  }
  for (int j = 0; j < res.y(); j++) {
    for (int i = 0; i < res.x(); i++) {
      std::cerr << sim[i][j][0] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << "===\n";
  return 0;
}