#include "levelset_fluid.h"
#include <iomanip>
#include <iostream>
#include <geometry/vector3.h>
#include <geometry/matrix.h>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <utils/file_writer.h>

int main(int argc, char const *argv[]) {
  LevelSetFluid3 sim;
  Ramuh::FileWriter writer;
  int resolution;
  std::stringstream folderName;
  //   writer.setDebug(true);
  if (argc < 2) {
    resolution = 32;
    folderName << "data" << resolution << '/';
  } else {
    // TODO: Read resolution
  }
  sim.setResolution(Ramuh::Vector3i(resolution));
  sim.setSize(Ramuh::Vector3d(1.0, 1.0, 1.0));

  std::cerr << "Initialized grid with size " << sim.domainSize()
            << ", resolution " << sim.resolution() << "  and h spacing "
            << sim.h() << std::endl;

  auto res = sim.resolution();
  auto h = sim.h();
  sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.7, 0.5), 0.15);
  // sim.addCubeSurface(Ramuh::Vector3d(0.45, 0.45, 0.45),
  //  Ramuh::Vector3d(0.75, 0.75, 0.75));
  sim.addCubeSurface(Ramuh::Vector3d(-15, -15, -15),
                     Ramuh::Vector3d(15, 0.3, 15));
  // sim.addCubeSurface(Ramuh::Vector3d(-5, -5, -5),
  //  Ramuh::Vector3d(0.2, 0.8, 0.2));

  //   sim.printLevelSetValue();
  //   std::cerr << std::endl;
  //   sim.redistance();
  //   sim.printLevelSetValue();
  sim.setVelocity();
  // writer.writeMeshModel(sim.marchingTetrahedra(), "data/model/0000.obj");
  sim.redistance();
  writer.writeMeshModel(sim.marchingTetrahedra(), "data/model/0000.obj");
  // sim.printVertexVelocity();
  // writer.writeLevelSet(sim, "data/0");
  for (int frame = 1; frame <= 100; frame++) {
    sim.checkCellMaterial();
    sim.addGravity();
    sim.boundaryVelocities();
    // if (frame >= 14)
    // sim.printFaceVelocity();
    sim.advectGridVelocity();
    // if (frame >= 14)
    // sim.printFaceVelocity();
    sim.solvePressure();
    // sim.printFaceVelocity();
    sim.extrapolateVelocity();
    // sim.printFaceVelocity();
    // sim.printLevelSetValue();
    sim.integrateLevelSet();
    // sim.printLevelSetValue();
    // std::cout << "plo" << std::endl;
    // if (!(frame % 5))
    sim.redistance();
    std::ostringstream filename;
    filename << "data/model/" << std::setw(4) << std::setfill('0') << frame
             << ".obj";
    writer.writeMeshModel(sim.marchingTetrahedra(), filename.str());
    // filename << "data128/" << frame;
    // writer.writeLevelSet(sim, std::string(filename.str()));
  }
  return 0;
}