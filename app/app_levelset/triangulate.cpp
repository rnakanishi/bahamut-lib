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
  sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.5, 0.5), 0.25);
  // sim.addCubeSurface(Ramuh::Vector3d(0.45, 0.45, 0.45),
  //  Ramuh::Vector3d(0.75, 0.75, 0.75));
  //   sim.addCubeSurface(Ramuh::Vector3d(-15, -15, -15),
  //  Ramuh::Vector3d(15, 0.3, 15));

  //   sim.redistance();

  writer.writeMeshModel(sim.marchingTetrahedra(), "obj/tetra.obj");

  writer.writeLevelSet(sim, "tetrahedrons");

  return 0;
}