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
#include <utils/timer.hpp>

int main(int argc, char const *argv[]) {
  LevelSetFluid3 sim;
  Ramuh::FileWriter writer;
  int resolution;
  std::stringstream folderName;
  //   writer.setDebug(true);
  if (argc < 2) {
    resolution = 20;
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
  sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.3, 0.5), 0.25);
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
  writer.writeMeshModel(sim.marchingTetrahedra(), "obj/0000.obj");
  // sim.printVertexVelocity();
  writer.writeLevelSet(sim, "data/0");

  Ramuh::Timer stopwatch;
  for (int frame = 1; frame <= 300; frame++) {
    Ramuh::TriangleMesh surface;
    try {
      sim.checkCellMaterial();
      sim.addGravity();
      sim.boundaryVelocities();

      stopwatch.reset();
      sim.macComarckVelocityAdvection();
      stopwatch.registerTime("Velocity advection");

      sim.solvePressure();
      stopwatch.registerTime("Pressure");

      sim.extrapolateVelocity();
      stopwatch.registerTime("Extraploate velocity");

      sim.macCormackAdvection();
      stopwatch.registerTime("Levelset advection");

      sim.redistance();
      stopwatch.registerTime("Redistance");

      surface = sim.marchingTetrahedra();
      stopwatch.registerTime("Marching Tetrahedra");
      stopwatch.evaluateComponentsTime();
    } catch (const char *error) {
      std::cerr << error << std::endl;
      return -1;
    }

    std::ostringstream filename, objname;
    objname << "obj/" << std::setw(4) << std::setfill('0') << frame << ".obj";
    writer.writeMeshModel(surface, objname.str());
    filename << "data/" << frame;
    writer.writeLevelSet(sim, std::string(filename.str()));
  }
  return 0;
}