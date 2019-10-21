#include <fluids/particle_levelset3.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utils/timer.hpp>
#include <utils/file_writer.h>
#include <pugixml.hpp>
#include <iterator>
#include <cstdio>

class PLevelSet3 : public Leviathan::ParticleLevelSet3 {
public:
  PLevelSet3() : ParticleLevelSet3() {}

  PLevelSet3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain)
      : ParticleLevelSet3(gridSize, domain) {
    _maxParticles = 80;
    lscount = 0;
    // _dt = 1 / 30.;
  }

  void initializeLevelSet(Eigen::Array3d center, double radius) {
    auto &function = getCellScalarData(_phiId);
#pragma omp parallel for
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getCellPosition(i);

      p -= center;
      // function[i] = -2;
      // if (sqrt(pow(p[0] - center[0], 2) + pow(p[1] - center[1], 2)) < radius)
      function[i] = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - radius;
      // function[i] = 5;
    }
  }

  void printParticles() {
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << baseFolder << "/" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);
    if (!file.is_open()) {
      std::cerr << "\033[1;31mError\033[0m Faile to open " << filename.str()
                << std::endl;
      return;
    }
    auto &vel = getParticleArrayData(_particleVelocityId);
    auto &signal = getParticleScalarData(_particleSignalId);
    auto &radius = getParticleScalarData(_particleRadiusId);

    for (size_t i = 0; i < getTotalParticleCount(); i++) {
      if (isActive(i)) {
        auto pos = getParticlePosition(i);
        file << pos[0] << " " << pos[1] << " " << pos[2] << " ";
        file << vel[i][0] << " " << vel[i][1] << " " << vel[i][2] << " ";
        file << signal[i] << " " << radius[i] << "\n";
      }
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

  void printLevelSet() {
    auto &phi = getCellScalarData(_phiId);
    std::ofstream file;
    std::stringstream filename;
    filename << baseFolder << "ls" << lscount++;
    file.open(filename.str().c_str(), std::ofstream::out);
    if (!file.is_open()) {
      std::cerr << "\033[1;31mError\033[0m Faile to open " << filename.str()
                << std::endl;
      return;
    }

    for (size_t i = 0; i < cellCount(); i++) {
      Eigen::Array3d pos = getCellPosition(i);
      file << pos[0] << " " << pos[1] << " " << pos[2] << " " << phi[i] << "\n";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

  void printGradients() {
    auto &gradient = getCellArrayData(_cellGradientId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << baseFolder << "/g" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);
    if (!file.is_open()) {
      std::cerr << "\033[1;31mError\033[0m Faile to open " << filename.str()
                << std::endl;
      return;
    }
    for (size_t i = 0; i < cellCount(); i++) {
      Eigen::Array3d pos = getCellPosition(i);
      file << pos[0] << " " << pos[1] << " " << pos[2] << " " << gradient[i][0]
           << " " << gradient[i][1] << " " << gradient[i][2] << "\n";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

  void printVelocities() {
    auto &velocity = getCellArrayData(_cellVelocityId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << baseFolder << "/v" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);
    if (!file.is_open()) {
      std::cerr << "\033[1;31mError\033[0m Faile to open " << filename.str()
                << std::endl;
      return;
    }
    for (size_t i = 0; i < cellCount(); i++) {
      Eigen::Array3d pos = getCellPosition(i);
      file << pos[0] << " " << pos[1] << " "
           << " " << pos[2] << " " << velocity[i][0] << " " << velocity[i][1]
           << " " << velocity[i][2] << "\n";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

  void defineVelocity(int i) {
    double T = 3.0;
    double time = (double)i * _originalDt / T;
    auto &u = getFaceScalarData(0, _cellVelocityId);
    auto &v = getFaceScalarData(1, _cellVelocityId);
    auto &w = getFaceScalarData(2, _cellVelocityId);

// u velocities
#pragma omp parallel
    {
#pragma omp for nowait
      for (size_t i = 0; i < faceCount(0); i++) {
        auto p = getFacePosition(0, i);
        u[i] = 2 * pow(sin(M_PI * p[0]), 2) * sin(2 * M_PI * p[1]) *
               sin(2 * M_PI * p[2]) * cos(M_PI * time);
        // u[i] *= ;

        // u[i] = 0;
      }

// v velocities
#pragma omp for nowait
      for (size_t j = 0; j < faceCount(1); j++) {
        auto p = getFacePosition(1, j);
        v[j] = -pow(sin(M_PI * p[1]), 2) * sin(2 * M_PI * p[0]) *
               sin(2 * M_PI * p[2]) * cos(M_PI * time);
        // v[j] *= ;
        // v[i] = -1;
      }

// w velocities
#pragma omp for nowait
      for (size_t k = 0; k < faceCount(2); k++) {
        auto p = getFacePosition(2, k);
        w[k] = -pow(sin(M_PI * p[2]), 2) * sin(2 * M_PI * p[0]) *
               sin(2 * M_PI * p[1]) * cos(M_PI * time);
        // w[k] *= ;
        // v[i] = -1;
      }
    }
  }
  void setBaseFolder(std::string folder) { baseFolder = folder; }

protected:
  int lscount;
  std::string baseFolder;
};

class Setup {
public:
  Setup() { _simulationType = -1; }

  Setup(std::string filename) : Setup() {
    auto parseResult = _document.load_file(filename.c_str());
    auto root = _document.children("simulation");

    _numberOfSimulations = std::distance(root.begin(), root.end());
    _currentSimulation = 0;

    getSimulationType();
  }

  void nextSimulation() {
    // TODO: ALL STRUCTURES: implement proper destructor for poiners
    _currentSimulation++;

    // Check simulation type
    getSimulationType();

    // Configure simulations parameters

    // run and save values
  }

  void prepareSimulation() {
    // Read parameters
    pugi::xml_node node;
    node = _document.child("simulation");
    for (int siblings = 0; siblings < _currentSimulation; siblings++)
      node = node.next_sibling();

    node = node.child("parameters");

    // Simulatino parameters
    int size;
    double domainMin, domainMax, dt = 1 / 60.;
    int maxParticles;

    sscanf(node.child("gridSize").child_value(), "%d", &size);
    sscanf(node.child("domain").child_value(), "%lf %lf", &domainMin,
           &domainMax);
    sscanf(node.child("frames").child_value(), "%d", &_frames);
    sscanf(node.child("dt").child_value(), "%lf", &dt);

    maxParticles = 80;
    if (!node.child("maxParticles").empty())
      sscanf(node.child("maxParticles").child_value(), "%d", &maxParticles);

    Eigen::Array3i gridSize(size);
    Ramuh::BoundingBox3 domain(domainMin, domainMax);

    _simulation = PLevelSet3(gridSize, domain);

    // Simulation output data
    std::string baseFolder(node.parent().child("baseFolder").child_value());
    _simulation.setBaseFolder(baseFolder);
    _logFilename =
        baseFolder + std::string(node.parent().child("logfile").child_value());

    // Confirming read data
    std::cout << "\n\n";
    std::cout << "**************************************\n";
    std::cout << "**************************************\n";
    std::cout << "Read the following data from xml file:\n";
    std::cout << "gridSize: " << size << std::endl;
    std::cout << "Frames: " << _frames << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "Base folder for results:  " << baseFolder << std::endl;

    // Initalizing domain
    Ramuh::Timer initTimer;

    _simulation.initializeLevelSet(Eigen::Array3d(.35, .35, .35), 0.15);
    _simulation.setMaxParticles(maxParticles);
    _simulation.setDt(dt);
    initTimer.registerTime("Initialization");
    _simulation.redistance();
    initTimer.registerTime("Redistance");

    _simulation.seedParticlesNearSurface();
    initTimer.registerTime("Particles seeding");
    _simulation.printParticles();
    _simulation.printLevelSet();
    initTimer.registerTime("printFile");
    initTimer.evaluateComponentsTime();
  }

  int getSimulationCount() { return _numberOfSimulations; }

  int getSimulationType() {
    pugi::xml_node node = _document.child("simulation");
    for (int siblings = 0; siblings < _currentSimulation; siblings++)
      node = node.next_sibling();

    node = node.child("parameters");
    node = node.child("type");

    _simulationTypeString = node.child_value();
    if (!_simulationTypeString.compare("particle_level_set"))
      _simulationType = 0;
    if (!_simulationTypeString.compare("weno_advection"))
      _simulationType = 1;
    if (!_simulationTypeString.compare("semi_lagrangean"))
      _simulationType = 2;

    if (_simulationType < 0)
      std::cerr << "\033[1;31m[ERROR]\033[0m"
                << " Parser failed on checking simulation type "
                << _document.name() << std::endl;

    std::cout << " Got " << _simulationTypeString
              << " type: " << _simulationType << "\n";

    return _simulationType;
  }

  void runAllSimulations() {
    for (size_t sim = 0; sim < _numberOfSimulations; sim++) {
      prepareSimulation();
      run();
      nextSimulation();
    }
  }

  void run() {

    pugi::xml_node node;
    node = _document.child("simulation");
    for (int siblings = 0; siblings < _currentSimulation; siblings++)
      node = node.next_sibling();
    std::string baseFolder(node.child("baseFolder").child_value());

    int lastRedistance = 0;
    int pReseed = 0;
    auto h = _simulation.getH();
    std::vector<int> particleHistory;
    std::vector<int> cellHistory;
    int advectionType = getSimulationType();

    Ramuh::Timer timer;
    for (int i = 0; i <= _frames; i++) {
      timer.clearAll();
      timer.reset();
      _simulation.defineVelocity(i);
      timer.registerTime("faceVelocity");
      _simulation.applyCfl();
      timer.registerTime("cfl");

      do {
        if (advectionType == 0 || advectionType == 2) {
          _simulation.findSurfaceCells(4.0 * h[0]);
          timer.registerTime("trackSurface");
        }

        if (advectionType == 2) {
          _simulation.computeCellVelocity();
          // _simulation.computeWenoGradient();
          // _simulation.advectUpwind();
          _simulation.advectSemiLagrangean();
        } else
          _simulation.advectWeno();
        timer.registerTime("cellAdvection");

        if (advectionType == 0) {
          _simulation.advectParticles((double)i / 3.0);
          timer.registerTime("particleAdvect");

          if (_simulation.correctLevelSetWithParticles()) {
            timer.registerTime("correction");

            lastRedistance = 0;
            _simulation.redistance();
            timer.registerTime("redistance");

            _simulation.correctLevelSetWithParticles();
            timer.registerTime("correction2");
          } else
            timer.registerTime("correction");
        }

        // if (lastRedistance >= 5) {
        lastRedistance = 0;
        _simulation.redistance();
        timer.registerTime("redistance");
        // }
        // lastRedistance++;
      } while (!_simulation.advanceTime());

      if (advectionType == 0) {
        // if (i % 5 == 0) {
        _simulation.reseedParticles();
        particleHistory.emplace_back(_simulation.getParticleCount());
        timer.registerTime("reseed");

        _simulation.adjustParticleRadius();
        timer.registerTime("radiusAdjust");
        std::cout << " Final count: " << _simulation.getParticleCount()
                  << " particles\n";
        // }
        // _simulation.printParticles();
      }

      _simulation.printLevelSet();
      timer.registerTime("print");

      timer.evaluateComponentsTime();
      if (i % 5 == 0)
        timer.evaluateComponentsAverageTime();

      timer.logToFile(_logFilename);
      Ramuh::FileWriter::writeArrayToFile(baseFolder + "particleCount.txt",
                                          particleHistory);

      cellHistory.emplace_back(_simulation.getSurfaceCellCount());
      Ramuh::FileWriter::writeArrayToFile(baseFolder + "cellCount.txt",
                                          cellHistory);
    }

    timer.evaluateComponentsAverageTime();
  }

private:
  int _simulationType;
  pugi::xml_document _document;
  int _numberOfSimulations, _currentSimulation;
  std::string _simulationTypeString;

  PLevelSet3 _simulation;
  int _frames;
  std::string _logFilename;
};

int main(int argc, char const *argv[]) {
  // PLevelSet3 system(Eigen::Array3i(100), Ramuh::BoundingBox3(0, 1));
  Setup configuration("./carbuncle/configs/pls.xml");
  configuration.runAllSimulations();

  return 0;
}
