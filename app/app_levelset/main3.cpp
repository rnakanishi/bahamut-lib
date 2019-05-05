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
#include <pugixml.hpp>

int main(int argc, char const *argv[]) {

  // Find a file to parse simulatoin parameters
  std::string xmlFilename;
  if (argc < 2) {
    std::cerr << "Expected a XML file\nUsing default config.xml instead\n";
    xmlFilename = std::string("assets/simulations/config.xml");
  } else {
    xmlFilename = std::string(argv[1]);
  }

  pugi::xml_document xmlDoc;
  if (!xmlDoc.load_file(xmlFilename.c_str())) {
    std::cerr << "Failed to load XML file: " << xmlFilename << std::endl;
    return -3;
  }

  LevelSetFluid3 sim;
  Ramuh::FileWriter writer;
  int resolution, nFrames;
  std::string objFolderName, dataFolderName;
  bool isPressureSecondOrder;

  // ================
  // Parsing values and assigning them
  pugi::xml_node node;
  node = xmlDoc.child("simulation");
  resolution = node.child("resolution").attribute("xyz").as_int();
  nFrames = node.child("maxFrames").attribute("value").as_int();
  dataFolderName = std::string(node.child("dataFolder").child_value());
  objFolderName = std::string(node.child("objFolder").child_value());
  int pressureOrder = node.child("pressureOrder").attribute("order").as_int();
  isPressureSecondOrder = (pressureOrder == 2) ? true : false;

  sim.setResolution(Ramuh::Vector3i(resolution));
  sim.setSize(Ramuh::Vector3d(1.0, 1.0, 1.0));
  sim.setPressureSecondOrder(isPressureSecondOrder);

  std::cerr << "Initialized grid with size " << sim.domainSize()
            << ", resolution " << sim.resolution() << "  and h spacing "
            << sim.h() << std::endl;
  std::cerr << "Pressure order at FS: " << ((isPressureSecondOrder) ? "2" : "1")
            << std::endl;

  auto res = sim.resolution();
  auto h = sim.h();

  sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.6, 0.5), 0.15);
  // sim.addSphereSurface(Ramuh::Vector3d(0.5, 0.55, 0.6), 0.25);
  // sim.addCubeSurface(Ramuh::Vector3d(-15, -15, -15),
  //  Ramuh::Vector3d(15, 0.4, 15));
  // sim.addCubeSurface(Ramuh::Vector3d(-5, -5, -5),
  //  Ramuh::Vector3d(0.2, 0.8, 0.2));

  sim.setVelocity();
  sim.redistance();
  std::ostringstream dataname, objname;

  dataname << dataFolderName << "/0";
  objname << objFolderName << "/0000.obj";

  writer.writeMeshModel(sim.marchingTetrahedra(), objname.str().c_str());
  writer.writeLevelSet(sim, dataname.str().c_str());

  Ramuh::Timer stopwatch;
  for (int frame = 1; frame <= nFrames; frame++) {
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

    objname.str(std::string());
    dataname.str(std::string());
    objname << objFolderName << "/" << std::setw(4) << std::setfill('0')
            << frame << ".obj";
    dataname << dataFolderName << "/" << frame;

    try {
      writer.writeMeshModel(surface, objname.str());
    } catch (const char *error) {
      std::cerr << "Failed to write obj: ";
      std::cerr << error << std::endl;
    }
    try {
      writer.writeLevelSet(sim, std::string(dataname.str()));
    } catch (const char *error) {
      std::cerr << "Failed to write levelset: " << error << std::endl;
    }
  }
  return 0;
}