#include <pugixml.hpp>
#include <Eigen/Dense>
#include <geometry/bounding_box.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utils/timer.hpp>
#include <surface/dual_cubes.h>
#include <utils/file_writer.h>

namespace Carbuncle {

class FileReader {
public:
  FileReader() { _simulationType = -1; }

  FileReader(std::string filename) : FileReader() {
    auto parseResult = _document.load_file(filename.c_str());
    auto root = _document.children("simulation");

    _numberOfSimulations = std::distance(root.begin(), root.end());
    _currentSimulation = 0;

    getSimulationType();
  }

  void nextSimulation() {
    _currentSimulation++;

    // Check simulation type
    getSimulationType();
  }

  void prepareSimulation(PLeveSet3 simulation) {
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
    double domainMin, domainMax, dt = 1 / 60.;
    int maxParticles;

    sscanf(node.child("gridSize").child_value(), "%d", &size);
    sscanf(node.child("domain").child_value(), "%lf %lf", &domainMin,
           &domainMax);
    sscanf(node.child("frames").child_value(), "%d", &_frames);
    sscanf(node.child("dt").child_value(), "%lf", &dt);

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

    maxParticles = 80;
    if (!node.child("maxParticles").empty()) {
      sscanf(node.child("maxParticles").child_value(), "%d", &maxParticles);
      configNode.append_child("maxParticles")
          .append_child(pugi::node_pcdata)
          .set_value(node.child("maxParticles").child_value());
    }

    Eigen::Array3i gridSize(size);
    Ramuh::BoundingBox3 domain(domainMin, domainMax);

    simulation = PLevelSet3(gridSize, domain);

    // Simulation output data
    node = node.parent();
    std::string baseFolder(node.child("baseFolder").child_value());
    simulation.setBaseFolder(baseFolder);
    std::string meshFolder(node.child("meshFolder").child_value());
    simulation.setMeshFolder(meshFolder);
    _logFilename =
        baseFolder + std::string(node.child("logfile").child_value());

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
    std::cout << "Base folder for results:  " << baseFolder << std::endl;
    if (!_doReseed)
      std::cout << "Particles won't be reseeded\n";

    // Initalizing domain
    Ramuh::Timer initTimer;

    simulation.initializeLevelSet(Eigen::Array3d(.35, .35, .35), 0.135);
    simulation.setMaxParticles(maxParticles);
    simulation.setDt(dt);
    initTimer.registerTime("Initialization");
    simulation.redistance();
    initTimer.registerTime("Redistance");

    simulation.findSurfaceCells(2.0);
    simulation.computeWenoGradient();
    initTimer.registerTime("Weno gradient");
    getSimulationType();
    if (_simulationType == 0 || _simulationType == 3) {
      simulation.seedParticlesNearSurface();
      initTimer.registerTime("Particles seeding");
      simulation.attractParticles();
      initTimer.registerTime("Particles attraction");
      simulation.printParticles();
    }
    // simulation.printLevelSet();
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
    auto h = simulation.getH();
    std::vector<int> particleHistory;
    std::vector<int> cellHistory;
    int advectionType = getSimulationType();

    Leviathan::DualCubes cubes(simulation);
    cubes.resetFileCounter();
    cubes.setFolder(simulation.getMeshFolder());
    cubes.computeIntersectionAndNormals();
    cubes.extractSurface();

    Ramuh::Timer timer;
    for (int i = 1; i <= _frames; i++) {
      timer.clearAll();
      timer.reset();

      simulation.defineVelocity(i);
      timer.registerTime("faceVelocity");
      simulation.applyCfl();
      timer.registerTime("cfl");

      do {
        if (advectionType != 3 && advectionType != 4) {
          simulation.findSurfaceCells(8.0 * h[0]);
          timer.registerTime("trackSurface");
        }

        if (advectionType == 2) {
          simulation.computeCellVelocity();
          // simulation.advectSemiLagrangeanThirdOrder();
          simulation.advectSemiLagrangean();
          timer.registerTime("cellAdvection");
        } else if (advectionType == 1 || advectionType == 0) {
          simulation.advectWeno();
          timer.registerTime("cellAdvection");
        } else if (advectionType == 4) {
          simulation.computeWenoGradient();
          simulation.advectCip();
          timer.registerTime("cellAdvection");
        }

        if (advectionType == 0) {
          simulation.advectParticles((double)i / 3.0);
          timer.registerTime("particleAdvect");
        }

        if (advectionType == 0) {
          if (simulation.correctLevelSetWithParticles()) {
            timer.registerTime("correction");

            lastRedistance = 0;
            simulation.redistance();
            timer.registerTime("redistance");

            simulation.correctLevelSetWithParticles();
            timer.registerTime("correction2");
          } else
            timer.registerTime("correction");
        }

        // lastRedistance++;
        if (advectionType != 3 && advectionType != 4) {
          // if (lastRedistance >= 5) {
          lastRedistance = 0;
          simulation.redistance();
          timer.registerTime("redistance");
          // }
        }
      } while (!simulation.advanceTime());

      if (advectionType == 0 && _doReseed) {
        if (i % 20 == 0) {
          simulation.reseedParticles();
          particleHistory.emplace_back(simulation.getParticleCount());
          timer.registerTime("reseed");

          simulation.adjustParticleRadius();
          timer.registerTime("radiusAdjust");
          std::cout << " Final count: " << simulation.getParticleCount()
                    << " particles\n";
        }
      }
      // if (advectionType == 0 || advectionType == 3)
      //   simulation.printParticles();

      // if (advectionType == 0 || advectionType == 3)
      //   simulation.printParticles();
      // if (advectionType != 3)
      // simulation.printLevelSet();
      // timer.registerTime("print");

      timer.logToFile(_logFilename);
      Ramuh::FileWriter::writeArrayToFile(baseFolder + "particleCount.txt",
                                          particleHistory);

      cellHistory.emplace_back(simulation.getSurfaceCellCount());
      Ramuh::FileWriter::writeArrayToFile(baseFolder + "cellCount.txt",
                                          cellHistory);

      cubes.swapLevelSet(simulation);
      cubes.resetFileCounter(i);
      cubes.setFolder(simulation.getMeshFolder());
      cubes.computeIntersectionAndNormals();
      cubes.extractSurface();
      timer.registerTime("extractSurface");
      timer.evaluateComponentsTime();

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

  bool _doReseed;
  int _frames;
  std::string _logFilename;
};

} // namespace Carbuncle