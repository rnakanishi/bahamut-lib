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
    resolution = 128;
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

  sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.5, 0.5), 0.15);
  // sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.55, 0.6), 0.25);
  sim.addCubeSurface(Ramuh::Vector3d(-15, -15, -15),
                     Ramuh::Vector3d(15, 0.2, 15));
  // sim.addCubeSurface(Ramuh::Vector3d(-5, -5, -5),
  //  Ramuh::Vector3d(0.2, 0.8, 0.2));

  sim.setVelocity();
  sim.redistance();
  writer.writeMeshModel(sim.marchingTetrahedra(), "obj/0000.obj");
  writer.writeLevelSet(sim, "data/0");

  sim.setPressureSecondOrder(false);
  Ramuh::Timer stopwatch;
  for (int frame = 1; frame <= 300; frame++) {
    std::cout << std::endl;
    Ramuh::TriangleMesh surface;
    try {
      stopwatch.reset();
      sim.checkCellMaterial();
      sim.addGravity();

      sim.extrapolateVelocity();
      stopwatch.registerTime("Gravity");

      stopwatch.registerTime("External forces");
      // sim.addExternalForce(Eigen::Vector3d(-1.5, -4.5, 0));
      sim.boundaryVelocities();

      sim.macComarckVelocityAdvection();
      // sim.advectGridVelocity();
      stopwatch.registerTime("Velocity advection");

      sim.writeVelocityField();
      std::cout << "Velocity advection\n";

      sim.extrapolateVelocity();
      stopwatch.registerTime("Extraploate velocity");

      sim.writeVelocityField();
      std::cout << "First extrapolation\n";

      sim.solvePressure();
      stopwatch.registerTime("Pressure");
      sim.boundaryVelocities();

      sim.extrapolateVelocity();
      stopwatch.registerTime("Extraploate velocity");

      sim.writeVelocityField();
      std::cout << "Second extrapolation\n";

      sim.macCormackAdvection();
      // sim.integrateLevelSet();
      stopwatch.registerTime("Levelset advection");

      if (frame % 5 == 0) {
        sim.redistance();
        stopwatch.registerTime("Redistance");
      }

      surface = sim.marchingTetrahedra();
      stopwatch.registerTime("Marching Tetrahedra");
      stopwatch.evaluateComponentsTime();
    } catch (const char *error) {
      std::cerr << error << std::endl;
      return -1;
    }

    std::ostringstream filename, objname;
    objname << "obj128/" << std::setw(4) << std::setfill('0') << frame
            << ".obj";
    filename << "data128/" << frame;
    try {
      writer.writeMeshModel(surface, objname.str());
    } catch (const char *error) {
      std::cerr << "Failed to write obj: ";
      std::cerr << error << std::endl;
    }
    try {
      writer.writeLevelSet(sim, std::string(filename.str()));
    } catch (const char *error) {
      std::cerr << "Failed to write levelset: " << error << std::endl;
    }
  }
  return 0;
}