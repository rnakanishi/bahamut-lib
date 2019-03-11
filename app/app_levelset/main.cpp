#include "levelset_fluid.h"
#include <iomanip>
#include <iostream>
#include <geometry/vector3.h>
#include <geometry/matrix.h>
#include <vector>
#include <cmath>
#include <utils/file_writer.h>

int main(void) {
  LevelSetFluid sim;
  Ramuh::FileWriter writer;
  // writer.setDebug(true);
  int resolution = 128;
  sim.setResolution(Ramuh::Vector3i(resolution, resolution, 1));
  sim.setSize(Ramuh::Vector3d(1.0, 1.0, 1.0 / resolution));

  std::cerr << "Initialized grid with size " << sim.domainSize()
            << ", resolution " << sim.resolution() << "  and h spacing "
            << sim.h() << std::endl;

  auto res = sim.resolution();
  auto h = sim.h();
  sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.7, 0), 0.1);
  sim.addCubeSurface(Ramuh::Vector3d(0, 0, 0),
                     Ramuh::Vector3d(1, 0.3, 2.0 / resolution));
  // sim.addCubeSurface(Ramuh::Vector3d(0, 0, 0),
  //  Ramuh::Vector3d(0.2, 0.8, 2.0 / resolution));

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

    std::ostringstream filename;
    // filename << "data/" << std::setw(4) << std::setfill('0') << frame;
    filename << "data/" << frame;
    writer.writeLevelSet(sim, std::string(filename.str()));
  }
  return 0;
}