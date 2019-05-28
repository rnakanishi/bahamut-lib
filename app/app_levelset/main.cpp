#include "../levelset_fluid.h"
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
    std::cerr << "Expected a XML file. Using default config.xml instead\n";
    xmlFilename = std::string("assets/simulations/config.xml");
  } else {
    xmlFilename = std::string(argv[1]);
  }

  LevelSetFluid2 sim;
  Ramuh::FileWriter writer;

  // =====================================================
  // =====================================================
  // Parsing values and assigning them
  sim.loadConfiguration(xmlFilename);

  auto res = sim.resolution();
  auto h = sim.h();
  // =====================================================
  // =====================================================

  sim.setVelocity();
  // sim.redistance();
  std::ostringstream dataname, objname;

  std::string &dataFolderName = sim.getDataFolderName();
  std::string &objFolderName = sim.getFolderName();

  dataname << sim.getDataFolderName() << "/0";
  objname << sim.getFolderName() << "/0000.obj";

  writer.writeMeshModel(sim.marchingTriangles(), objname.str().c_str());
  writer.writeLevelSet(sim, dataname.str().c_str());

  int nFrames = sim.getFramesNumber();
  int velocityAdvectionOrder = sim.getVelocityOrder();
  int levelsetAdvectionOrder = sim.getLevelsetOrder();

  for (int frame = 1; frame <= nFrames; frame++) {
    std::cout << std::endl;
    Ramuh::TriangleMesh surface;
    sim.run();

    objname.str(std::string());
    dataname.str(std::string());
    objname << objFolderName << "/" << std::setw(4) << std::setfill('0')
            << frame << ".obj";
    dataname << dataFolderName << "/" << frame;

    try {
      surface = sim.marchingTriangles();
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