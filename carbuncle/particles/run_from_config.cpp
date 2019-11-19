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
#include <surface/dual_cubes.h>
#include <omp.h>

class PLevelSet3 : public Leviathan::ParticleLevelSet3 {
public:
  PLevelSet3() : ParticleLevelSet3() {}

  PLevelSet3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain)
      : ParticleLevelSet3(gridSize, domain) {
    _maxParticles = 80;
    lscount = 0;
    pfileCount = 0;

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

  void initializeGradient(Eigen::Array3d center, double radius) {
    auto &function = getCellArrayData(_cellGradientId);
#pragma omp parallel for
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getCellPosition(i);
      p -= center;
      for (size_t dim = 0; dim < 3; dim++)
        function[i][dim] = 2 * p[dim];
      // function[i] = function[i].matrix().normalized().array();
      // int a = 0;
    }
  }

  void printParticles() {
    std::ofstream fileneg;
    std::ofstream filepos;
    std::stringstream negFile;
    std::stringstream posFile;
    negFile << baseFolder << "pneg" << pfileCount << ".obj";
    posFile << baseFolder << "ppos" << pfileCount++ << ".obj";
    fileneg.open(negFile.str().c_str(), std::ofstream::out);
    filepos.open(posFile.str().c_str(), std::ofstream::out);
    if (!fileneg.is_open() || !filepos.is_open()) {
      std::cerr << "\033[1;31mError\033[0m Faile to open " << negFile.str()
                << std::endl;
      return;
    }
    auto &vel = getParticleArrayData(_particleVelocityId);
    auto &signal = getParticleScalarData(_particleSignalId);
    auto &radius = getParticleScalarData(_particleRadiusId);

    for (size_t i = 0; i < getTotalParticleCount(); i++) {
      if (isActive(i)) {
        auto pos = getParticlePosition(i);
        if (signal[i] < 0) {
          fileneg << "v ";
          fileneg << pos[0] << " " << pos[1] << " " << pos[2] << " ";
          // fileneg << vel[i][0] << " " << vel[i][1] << " " << vel[i][2] << "
          // "; fileneg << signal[i]; " " << radius[i];
          fileneg << std::endl;
        } else {
          filepos << "v ";
          filepos << pos[0] << " " << pos[1] << " " << pos[2] << " ";
          // filepos << vel[i][0] << " " << vel[i][1] << " " << vel[i][2] << "
          // "; filepos << signal[i]; " " << radius[i];
          filepos << std::endl;
        }
      }
    }
    std::cout << "File written: " << negFile.str() << std::endl;
    fileneg.close();
    filepos.close();
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
        // u[i] = p[1];
      }

// v velocities
#pragma omp for nowait
      for (size_t j = 0; j < faceCount(1); j++) {
        auto p = getFacePosition(1, j);
        v[j] = -pow(sin(M_PI * p[1]), 2) * sin(2 * M_PI * p[0]) *
               sin(2 * M_PI * p[2]) * cos(M_PI * time);
        // v[j] = -p[0];
      }

// w velocities
#pragma omp for nowait
      for (size_t k = 0; k < faceCount(2); k++) {
        auto p = getFacePosition(2, k);
        w[k] = -pow(sin(M_PI * p[2]), 2) * sin(2 * M_PI * p[0]) *
               sin(2 * M_PI * p[1]) * cos(M_PI * time);
        // w[k] = 0;
      }
    }
  }

  void setBaseFolder(std::string folder) { baseFolder = folder; }

  void setMeshFolder(std::string folder) { meshFolder = folder; }

  std::string getMeshFolder() { return meshFolder; }

protected:
  int lscount;
  int pfileCount;
  std::string baseFolder, meshFolder;
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
    pugi::xml_document configDoc;

    node = _document.child("simulation");
    for (int siblings = 0; siblings < _currentSimulation; siblings++)
      node = node.next_sibling();

    auto configNode = configDoc.append_child("simulation");
    configNode = configNode.append_child("parameter");
    node = node.child("parameters");

    // Simulatino parameters
    int size;
    double domainMin, domainMax, dt = 1 / 60., cflCond = -1;
    int maxParticles;

    sscanf(node.child("gridSize").child_value(), "%d", &size);
    sscanf(node.child("domain").child_value(), "%lf %lf", &domainMin,
           &domainMax);
    sscanf(node.child("frames").child_value(), "%d", &_frames);
    sscanf(node.child("dt").child_value(), "%lf", &dt);
    sscanf(node.child("dt").child_value(), "%lf", &dt);
    sscanf(node.child("cfl").child_value(), "%lf", &cflCond);

    configNode.append_child("gridSize")
        .append_child(pugi::node_pcdata)
        .set_value(node.child("gridSize").child_value());
    configNode.append_child("domain")
        .append_child(pugi::node_pcdata)
        .set_value(node.child("domain").child_value());
    configNode.append_child("frames")
        .append_child(pugi::node_pcdata)
        .set_value(node.child("frames").child_value());
    configNode.append_child("dt")
        .append_child(pugi::node_pcdata)
        .set_value(node.child("dt").child_value());
    if (cflCond > 0)
      configNode.append_child("cfl")
          .append_child(pugi::node_pcdata)
          .set_value(node.child("cfl").child_value());

    maxParticles = 80;
    if (!node.child("maxParticles").empty()) {
      sscanf(node.child("maxParticles").child_value(), "%d", &maxParticles);
      configNode.append_child("maxParticles")
          .append_child(pugi::node_pcdata)
          .set_value(node.child("maxParticles").child_value());
    }

    Eigen::Array3i gridSize(size);
    Ramuh::BoundingBox3 domain(domainMin, domainMax);

    _simulation = PLevelSet3(gridSize, domain);

    // Simulation output data
    node = node.parent();
    std::string baseFolder(node.child("baseFolder").child_value());
    _simulation.setBaseFolder(baseFolder);
    std::string meshFolder(node.child("meshFolder").child_value());
    _simulation.setMeshFolder(meshFolder);
    _logFilename =
        baseFolder + std::string(node.child("logfile").child_value());
    _statisticsFile =
        baseFolder + std::string(node.child("statisticsFile").child_value());

    configNode = configNode.parent();
    configNode.append_child("baseFolder")
        .append_child(pugi::node_pcdata)
        .set_value(node.child("baseFolder").child_value());
    configNode.append_child("meshFolder")
        .append_child(pugi::node_pcdata)
        .set_value(node.child("meshFolder").child_value());
    configNode.append_child("logfile")
        .append_child(pugi::node_pcdata)
        .set_value(node.child("logfile").child_value());

    _doReseed = true;
    node = node.child("meta");
    if (node) {
      std::string reseed(node.attribute("reseed").value());
      if (!reseed.compare("false")) {
        _doReseed = false;
        configNode.append_child("meta").append_attribute("reseed").set_value(
            node.attribute("reseed").value());
      }
    }
    std::ofstream xmlConfigFile;
    std::string xmlFilename(baseFolder + "config.xml");
    xmlConfigFile.open(xmlFilename);
    if (!xmlConfigFile.is_open())
      std::cerr << "Error opening " << xmlFilename << std::endl;
    else {
      configDoc.print(xmlConfigFile);
      xmlConfigFile.close();
    }

    // Confirming read data
    std::cout << "\n\n";
    std::cout << "**************************************\n";
    std::cout << "**************************************\n";
    std::cout << "Read the following data from xml file:\n";
    std::cout << "gridSize: " << size << std::endl;
    std::cout << "Frames: " << _frames << std::endl;
    std::cout << "dt: " << dt << std::endl;
    if (cflCond > 0)
      std::cout << "CFL condition: " << cflCond << std::endl;
    std::cout << "Base folder for results:  " << baseFolder << std::endl;
    if (!_doReseed)
      std::cout << "Particles won't be reseeded\n";

    // Initalizing domain
    Ramuh::Timer initTimer;

    _simulation.initializeLevelSet(Eigen::Array3d(.35, .35, .35), 0.15);
    _simulation.setMaxParticles(maxParticles);
    _simulation.setDt(dt);
    if (cflCond > 0)
      _simulation.setCflCondition(cflCond);
    initTimer.registerTime("Initialization");

    _simulation.setMaxRedistanceIterations(100);
    _simulation.redistance();
    _simulation.setMaxRedistanceIterations(15);
    initTimer.registerTime("Redistance");

    _simulation.findSurfaceCells(2.0);
    // _simulation.computeWenoGradient();
    _simulation.initializeGradient(Eigen::Array3d(.35, .35, .35), 0.15);
    initTimer.registerTime("Weno gradient");
    getSimulationType();
    if (_simulationType == 0 || _simulationType == 3) {
      _simulation.seedParticlesNearSurface();
      initTimer.registerTime("Particles seeding");
      _simulation.attractParticles();
      initTimer.registerTime("Particles attraction");
      _simulation.printParticles();
    }
    // _simulation.printLevelSet();
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
    else if (!_simulationTypeString.compare("weno_advection"))
      _simulationType = 1;
    else if (!_simulationTypeString.compare("semi_lagrangean"))
      _simulationType = 2;
    else if (!_simulationTypeString.compare("particles"))
      _simulationType = 3;
    else if (!_simulationTypeString.compare("cip_advection"))
      _simulationType = 4;

    if (_simulationType < 0)
      std::cerr << "\033[1;31m[ERROR]\033[0m"
                << " Parser failed on checking simulation type "
                << _document.name() << std::endl;

    std::cout << " Got " << _simulationTypeString
              << " type: " << _simulationType << "\n";

    return _simulationType;
  }

  void extractStatistics(Ramuh::Statistics &statistics, Ramuh::Timer &timer) {
    auto h = _simulation.getH();
    auto surface = _simulation.findSurfaceCells(9.0 * h[0]);
    statistics.registerComponent("bandCells", (int)surface.size());
    statistics.registerComponent("fluidCells", _simulation.fluidCellCount());

    std::vector<std::string> fields;
    fields.emplace_back("faceVelocity");
    fields.emplace_back("cfl");
    fields.emplace_back("cellAdvection");
    fields.emplace_back("extractSurface");
    fields.emplace_back("total");
    fields.emplace_back("trackSurface");
    fields.emplace_back("redistace");
    switch (_simulationType) {
    case 0:
      fields.emplace_back("particleAdvect");
      fields.emplace_back("correction");
      fields.emplace_back("correction2");
      fields.emplace_back("redistance");
      fields.emplace_back("radiusAdjust");
      fields.emplace_back("reseed");
      break;
    case 1:
      break;
    case 2:
      break;
    default:
      break;
    }

    for (auto &field : fields)
      statistics.registerComponent(field, timer.getComponentTime(field));
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

    Leviathan::DualCubes cubes(_simulation);
    cubes.resetFileCounter();
    cubes.setFolder(_simulation.getMeshFolder());
    cubes.computeIntersectionAndNormals();
    cubes.extractSurface();

    Ramuh::Statistics statistics;
    Ramuh::Timer timer;
    for (int i = 1; i <= _frames; i++) {
      timer.clearAll();
      timer.reset();

      _simulation.defineVelocity(i);
      timer.registerTime("faceVelocity");
      int cflStep = _simulation.applyCfl();
      statistics.registerComponent("cflSteps", cflStep);
      timer.registerTime("cfl");

      do {
        if (advectionType != 3 && advectionType != 4) {
          auto surface = _simulation.findSurfaceCells(9.0 * h[0]);
          timer.registerTime("trackSurface");
          // _simulation.computeCellsGradient();
          // timer.registerTime("cellGradient");
        }

        if (advectionType == 2) {
          _simulation.computeCellVelocity();
          // _simulation.advectSemiLagrangeanThirdOrder();
          _simulation.advectSemiLagrangean();
          timer.registerTime("cellAdvection");
        } else if (advectionType == 1 || advectionType == 0) {
          _simulation.advectWeno();
          timer.registerTime("cellAdvection");
        } else if (advectionType == 4) {
          // _simulation.findSurfaceCells(2.0);
          // timer.registerTime("trackSurface");
          // _simulation.computeWenoGradient();
          // timer.registerTime("cellGradient");
          _simulation.advectCip();
          timer.registerTime("cellAdvection");
          // _simulation.computeWenoGradient();
        }

        if (advectionType == 0) {
          _simulation.advectParticles((double)i / 3.0);
          timer.registerTime("particleAdvect");
        }

        if (advectionType == 0) {
          if (_simulation.correctLevelSetWithParticles()) {
            timer.registerTime("correction");

            lastRedistance = 0;
            _simulation.redistance(true);
            timer.registerTime("redistance");

            _simulation.correctLevelSetWithParticles();
            timer.registerTime("correction2");

            _simulation.adjustParticleRadius();
            timer.registerTime("radiusAdjust");
          } else
            timer.registerTime("correction");
        }

        if (advectionType != 3 && advectionType != 5) {
          lastRedistance++;
          if (lastRedistance >= 6) {
            lastRedistance = 0;
            _simulation.redistance();
            timer.registerTime("redistance");
          }
        }
      } while (!_simulation.advanceTime());

      if (advectionType == 0 && _doReseed) {
        if (i % 25 == 0) {
          _simulation.reseedParticles();
          particleHistory.emplace_back(_simulation.getParticleCount());
          timer.registerTime("reseed");

          _simulation.adjustParticleRadius();
          timer.registerTime("radiusAdjust");
          std::cout << " Final count: " << _simulation.getParticleCount()
                    << " particles\n";
        }
      }
      // if (advectionType == 0 || advectionType == 3)
      //   _simulation.printParticles();

      // if (advectionType == 0 || advectionType == 3)
      //   _simulation.printParticles();
      // if (advectionType != 3)
      // _simulation.printLevelSet();
      // timer.registerTime("print");

      timer.logToFile(_logFilename);
      Ramuh::FileWriter::writeArrayToFile(baseFolder + "particleCount.txt",
                                          particleHistory);

      cellHistory.emplace_back(_simulation.getSurfaceCellCount());
      Ramuh::FileWriter::writeArrayToFile(baseFolder + "cellCount.txt",
                                          cellHistory);

      cubes.swapLevelSet(_simulation);
      cubes.resetFileCounter(i);
      cubes.setFolder(_simulation.getMeshFolder());
      // cubes.computeWenoGradient();
      cubes.computeIntersectionAndNormals();
      cubes.extractSurface();
      timer.registerTime("extractSurface");
      timer.evaluateComponentsTime();

      extractStatistics(statistics, timer);
      statistics.writeToFile(_statisticsFile);

      if (i % 5 == 0)
        timer.evaluateComponentsAverageTime();
      std::cout << "Finished timestep " << i << std::endl << std::endl;
    }

    timer.evaluateComponentsAverageTime();
  }

private:
  int _simulationType;
  pugi::xml_document _document;
  int _numberOfSimulations, _currentSimulation;
  std::string _simulationTypeString;

  PLevelSet3 _simulation;
  bool _doReseed;
  int _frames;
  std::string _logFilename, _statisticsFile;
};

int main(int argc, char const *argv[]) {
  omp_set_num_threads(14);
  // PLevelSet3 system(Eigen::Array3i(100), Ramuh::BoundingBox3(0, 1));
  Setup configuration("./carbuncle/configs/pls.xml");
  configuration.runAllSimulations();

  return 0;
}
