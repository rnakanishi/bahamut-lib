#include "levelset_fluid.h"
#include <iomanip>
#include <iostream>
#include <geometry/vector3.h>
#include <geometry/matrix.h>
#include <vector>
#include <cmath>
#include <utils/file_writer.h>

int main(void) {
  LevelSetFluid2 sim;
  Ramuh::FileWriter writer;
  // writer.setDebug(true);
  int resolution = 32;
  sim.setResolution(Ramuh::Vector2i(resolution, resolution));
  sim.setSize(Ramuh::Vector2d(1.0, 1.0));

  std::cerr << "Initialized grid with size " << sim.domainSize()
            << ", resolution " << sim.resolution() << "  and h spacing "
            << sim.h() << std::endl;

  auto res = sim.resolution();
  auto h = sim.h();
  sim.addSphereSurface(Ramuh::Vector2d(0.5, 0.8), 0.15);
  sim.addCubeSurface(Ramuh::Vector2d(-5, -5), Ramuh::Vector2d(5, 0.37));
  // sim.addCubeSurface(Ramuh::Vector2d(0, 0, 0),
  //  Ramuh::Vector2d(0.2, 0.8, 2.0 / resolution));

  // sim.redistance();
  sim.setVelocity();
  // sim.printVertexVelocity();
  writer.writeLevelSet(sim, "data/0");
  for (int frame = 1; frame < 150; frame++) {
    sim.checkCellMaterial();
    sim.addGravity();
    sim.advectGridVelocity();
    sim.boundaryVelocities();
    sim.solvePressure();
    // sim.printFaceVelocity();
    sim.extrapolateVelocity();
    // sim.printFaceVelocity();
    sim.integrateLevelSet();
    // if (!(frame % 5))
    //   sim.redistance();
    std::ostringstream filename;
    // filename << "data/" << std::setw(4) << std::setfill('0') << frame;
    filename << "data/" << frame;
    writer.writeLevelSet(sim, std::string(filename.str()));
  }
  return 0;
}