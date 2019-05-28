#include "levelset_fluid.h"
#include <fstream>
#include <geometry/vector3.h>
#include <pugixml.hpp>
#include <utils/timer.hpp>

LevelSetFluid2::LevelSetFluid2() : Ramuh::LevelSet2() {}

void LevelSetFluid2::run() {
  Ramuh::Timer stopwatch;
  // Correct cfl timestep
  cfl();

  // Run through timesteps given by cfl condition
  do {
    stopwatch.reset();
    stopwatch.clearAll();
    try {
      checkCellMaterial();
      addGravity();
      // addExternalForce(Eigen::Vector3d(0, -2.5, 0));
      stopwatch.registerTime("Gravity");

      extrapolateVelocity();
      boundaryVelocities();

      if (_velocityAdvectionOrder == 1)
        advectGridVelocity();
      else if (_velocityAdvectionOrder == 2)
        macCormackVelocityAdvection();
      stopwatch.registerTime("Velocity advection");

      writeVelocityField();
      std::cout << "Velocity advection\n";

      extrapolateVelocity();
      stopwatch.registerTime("Extrapolate velocity");
      boundaryVelocities();

      solvePressure();
      stopwatch.registerTime("Pressure");
      boundaryVelocities();
      writeVelocityField();

      extrapolateVelocity();
      boundaryVelocities();
      stopwatch.registerTime("Extrapolate velocity");

      if (_levelsetAdvectionOrder == 1)
        integrateLevelSet();
      else if (_levelsetAdvectionOrder == 2)
        macCormackAdvection();
      stopwatch.registerTime("Levelset advection");

      // if (frame % 5 == 0) {
      redistance();
      stopwatch.registerTime("Redistance");
      // }

      stopwatch.evaluateComponentsTime();
      // writeFaceVelocity("results/lastVelocity");
      // writeLevelSetValue("results/lastLevelset");
    } catch (const char *error) {
      std::cerr << error << std::endl;
      throw("Simulation error\n");
    }
  } while (!advanceTime());
}

void LevelSetFluid2::writeVelocityField() {
  std::ofstream file;
  file.open("results/velocityField");
  if (file.is_open()) {
    for (int id = 0; id < cellCount(); id++) {
      Eigen::Array2i ij = idToij(id);
      int i, j;
      i = ij[0];
      j = ij[1];

      file << (_u[i][j] + _u[i + 1][j]) / 2 << " ";
      file << (_v[i][j] + _v[i][j + 1]) / 2 << " ";
    }
  }
}

void LevelSetFluid2::loadConfiguration(std::string filename) {

  pugi::xml_document xmlDoc;
  pugi::xml_parse_result result = xmlDoc.load_file(filename.c_str());
  if (!result) {
    std::cerr << "Failed to load XML file: " << filename << std::endl;
    std::cerr << result.description() << std::endl;

    return;
  }

  pugi::xml_node node;
  node = xmlDoc.child("simulation");
  _resolution = node.child("resolution").attribute("xyz").as_int();
  _nFrames = node.child("maxFrames").attribute("value").as_int();
  _dataFolderName = std::string(node.child("dataFolder").child_value());
  _objFolderName = std::string(node.child("objFolder").child_value());
  int pressureOrder = node.child("pressure").attribute("order").as_int();
  _velocityAdvectionOrder =
      node.child("velocityAdvection").attribute("order").as_int();
  _levelsetAdvectionOrder =
      node.child("levelsetAdvection").attribute("order").as_int();
  _isPressure2nd = (pressureOrder == 2) ? true : false;

  setResolution(Eigen::Array2i(_resolution, _resolution));
  setSize(Eigen::Array2d(1.0, 1.0));
  setPressureOrder(_isPressure2nd);

  std::cerr << "Initialized grid with size " << domainSize() << ", resolution "
            << resolution() << "  and h spacing " << h() << std::endl;
  std::cerr << "Pressure order at FS: " << ((_isPressure2nd) ? "2" : "1")
            << std::endl;

  pugi::xml_node domain = node.child("domain");
  for (pugi::xml_node sphere : domain.children("sphere")) {
    Eigen::Array2d position;
    double radius;
    std::stringstream line;
    line << sphere.child("radius").child_value() << " ";
    line << sphere.child("center").child_value();

    line >> radius;
    line >> position[0] >> position[1];
    std::cerr << "Read sphere at ";
    std::cerr << position.transpose();
    std::cerr << " with radius " << radius << std::endl;
    addSphereSurface(position, radius);
  }
  for (pugi::xml_node cube : domain.children("cube")) {
    Eigen::Array2d lower, upper;
    double skip;
    std::stringstream line;
    line << cube.child("lower").child_value() << " ";
    line << cube.child("upper").child_value();

    line >> lower[0] >> lower[1] >> skip;
    line >> upper[0] >> upper[1];
    std::cerr << "Read cube from " << lower.transpose();
    std::cerr << " to " << upper.transpose() << std::endl;
    addCubeSurface(lower, upper);
  }
}

std::string &LevelSetFluid2::getFolderName() { return _objFolderName; }

std::string &LevelSetFluid2::getDataFolderName() { return _dataFolderName; }

int LevelSetFluid2::getVelocityOrder() { return (_velocityAdvectionOrder); }

int LevelSetFluid2::getLevelsetOrder() { return _levelsetAdvectionOrder; }

int LevelSetFluid2::getFramesNumber() { return _nFrames; }
