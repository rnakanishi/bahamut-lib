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
    std::cerr << "Expected a XML file. Using default config.xml instead\n";
    xmlFilename = std::string("assets/simulations/config.xml");
  } else {
    xmlFilename = std::string(argv[1]);
  }

  pugi::xml_document xmlDoc;
  pugi::xml_parse_result result = xmlDoc.load_file(xmlFilename.c_str());
  if (!result) {
    std::cerr << "Failed to load XML file: " << xmlFilename << std::endl;
    std::cerr << result.description() << std::endl;

    return -3;
  }

  LevelSetFluid3 sim;
  Ramuh::FileWriter writer;
  int resolution, nFrames;
  std::string objFolderName, dataFolderName;
  bool isPressureSecondOrder;

  // =====================================================
  // =====================================================
  // TODO: implement proper Parser class
  // Parsing values and assigning them
  pugi::xml_node node;
  node = xmlDoc.child("simulation");
  resolution = node.child("resolution").attribute("xyz").as_int();
  nFrames = node.child("maxFrames").attribute("value").as_int();
  dataFolderName = std::string(node.child("dataFolder").child_value());
  objFolderName = std::string(node.child("objFolder").child_value());
  int pressureOrder = node.child("pressure").attribute("order").as_int();
  int advectionOrder = node.child("advection").attribute("order").as_int();
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

  pugi::xml_node domain = node.child("domain");
  for (pugi::xml_node sphere : domain.children("sphere")) {
    Eigen::Array3d position;
    double radius;
    std::stringstream line;
    line << sphere.child("radius").child_value() << " ";
    line << sphere.child("center").child_value();

    line >> radius;
    line >> position[0] >> position[1] >> position[2];
    std::cerr << "Read sphere at ";
    std::cerr << position.transpose();
    std::cerr << " with radius " << radius << std::endl;
    sim.addSphereSurface(position, radius);
  }
  for (pugi::xml_node sphere : domain.children("cube")) {
    Eigen::Array3d lower, upper;
    std::stringstream line;
    line << sphere.child("lower").child_value() << " ";
    line << sphere.child("upper").child_value();

    line >> lower[0] >> lower[1] >> lower[2];
    line >> upper[0] >> upper[1] >> upper[2];
    std::cerr << "Read cube from " << lower.transpose();
    std::cerr << " to " << upper.transpose() << std::endl;
    sim.addCubeSurface(lower, upper);
  }
  // =====================================================
  // =====================================================

  sim.setVelocity();
  // sim.redistance();
  std::ostringstream dataname, objname;

  dataname << dataFolderName << "/0";
  objname << objFolderName << "/0000.obj";

  writer.writeMeshModel(sim.marchingTetrahedra(), objname.str().c_str());
  writer.writeLevelSet(sim, dataname.str().c_str());

  Ramuh::Timer stopwatch;
  for (int frame = 1; frame <= nFrames; frame++) {
    std::cout << std::endl;
    Ramuh::TriangleMesh surface;
    sim.cfl();

    stopwatch.reset();
    stopwatch.clearAll();
    do {
      try {
        sim.checkCellMaterial();
        sim.addGravity();
        // sim.addExternalForce(Eigen::Vector3d(-1.5, -4.5, 0));
        stopwatch.registerTime("Gravity");

        sim.extrapolateVelocity();
        sim.boundaryVelocities();

        if (advectionOrder == 1)
          sim.advectGridVelocity();
        else if (advectionOrder == 2)
          sim.macComarckVelocityAdvection();
        stopwatch.registerTime("Velocity advection");

        sim.writeVelocityField();
        std::cout << "Velocity advection\n";

        sim.extrapolateVelocity();
        stopwatch.registerTime("Extrapolate velocity");

        sim.solvePressure();
        stopwatch.registerTime("Pressure");
        // sim.boundaryVelocities();

        sim.extrapolateVelocity();
        stopwatch.registerTime("Extrapolate velocity");

        if (advectionOrder == 1)
          sim.integrateLevelSet();
        else if (advectionOrder == 2)
          sim.macCormackAdvection();
        stopwatch.registerTime("Levelset advection");

        // if (frame % 5 == 0) {
        // sim.redistance();
        // stopwatch.registerTime("Redistance");
        // }

        stopwatch.evaluateComponentsTime();
        sim.writeFaceVelocity("results/lastVelocity");
        sim.writeLevelSetValue("results/lastLevelset");
      } catch (const char *error) {
        std::cerr << error << std::endl;
        return -1;
      }
    } while (!sim.advanceTime());

    objname.str(std::string());
    dataname.str(std::string());
    objname << objFolderName << "/" << std::setw(4) << std::setfill('0')
            << frame << ".obj";
    dataname << dataFolderName << "/" << frame;

    try {
      surface = sim.marchingTetrahedra();
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