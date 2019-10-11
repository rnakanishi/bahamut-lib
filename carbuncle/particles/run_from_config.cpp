#include <fluids/particle_levelset3.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utils/timer.hpp>
#include <utils/file_writer.h>
#include <pugixml.hpp>

class PLevelSet3 : public Leviathan::ParticleLevelSet3 {
public:
  PLevelSet3() : ParticleLevelSet3() {}

  PLevelSet3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain)
      : ParticleLevelSet3(gridSize, domain) {
    _maxParticles = 80;
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
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << baseFolder << "/ls" << count++;
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
  std::string baseFolder;
};

class Setup {
public:
  Setup() { _simulationType = -1; }

  Setup(std::string filename) : Setup() {
    pugi::xml_document xmlDoc;
    _parseResult = xmlDoc.load_file(filename.c_str());
    pugi::xml_node node;
    node = xmlDoc.child("simulation");
    node = node.child("type");

    std::string type(node.child_value());
    if (!type.compare("particle_level_set"))
      _simulationType = 0;
    if (!type.compare("weno_advection"))
      _simulationType = 1;
    if (!type.compare("semi_lagrangean"))
      _simulationType = 2;

    if (_simulationType < 0)
      std::cerr << "\033[1;31m[ERROR]\033[0m"
                << " Parser failed on checking simulation type " << filename
                << std::endl;

    std::cerr << " Got " << type << " type " << _simulationType << "\n";
  }

  int getSimulationType() {
    // 0 Particle levelset
    // 1 weno Advection
    // 2 Semi lagrangean advection
    return _simulationType;
  }

  std::string getSimulationName() {
    if (_simulationType == 0)
      return "Particle Level Set";
    if (_simulationType == 1)
      return "WENO advection scheme";
    if (_simulationType == 2)
      return "Semi Lagrangean scheme";
    return "";
  }

  void configureSimulation(PLevelSet3 &simulation) {}

private:
  int _simulationType;
  pugi::xml_parse_result _parseResult;
};

int main(int argc, char const *argv[]) {
  PLevelSet3 system(Eigen::Array3i(100), Ramuh::BoundingBox3(0, 1));
  Setup configuration("./carbuncle/configs/pls.xml");
  system.setBaseFolder("results/particles/3d/pls_deform");

  Ramuh::Timer initTimer;

  system.initializeLevelSet(Eigen::Array3d(.35, .35, .35), 0.15);
  initTimer.registerTime("Initialization");
  system.redistance();
  initTimer.registerTime("Redistance");

  system.seedParticlesNearSurface();
  initTimer.registerTime("Particles seeding");
  system.printParticles();
  system.printLevelSet();
  initTimer.registerTime("printFile");
  initTimer.evaluateComponentsTime();

  int lastRedistance = 0;
  int pReseed = 0;
  auto h = system.getH();
  std::vector<int> particleHistory;
  Ramuh::Timer timer;
  std::cerr << "Starting simulation: " << configuration.getSimulationName()
            << "\n";
  int advectionType = configuration.getSimulationType();

  for (int i = 0; i <= 180; i++) {
    timer.clearAll();
    timer.reset();
    system.defineVelocity(i);
    timer.registerTime("faceVelocity");
    system.applyCfl();
    timer.registerTime("cfl");

    do {
      if (advectionType == 0 || advectionType == 1) {
        system.findSurfaceCells(4.0 * h[0]);
        timer.registerTime("trackSurface");
      }

      if (advectionType == 2) {
        system.computeCellVelocity();
        // system.computeWenoGradient();
        // system.advectUpwind();
        system.advectSemiLagrangean();
      } else
        system.advectWeno();
      timer.registerTime("cellAdvection");

      if (advectionType == 0) {
        system.advectParticles((double)i / 3.0);
        timer.registerTime("particleAdvect");
        particleHistory.emplace_back(system.getParticleCount());

        if (system.correctLevelSetWithParticles()) {
          timer.registerTime("correction");

          lastRedistance = 0;
          system.redistance();
          timer.registerTime("redistance");

          system.correctLevelSetWithParticles();
          timer.registerTime("correction2");
        } else
          timer.registerTime("correction");
      }

      // if (lastRedistance >= 5) {
      lastRedistance = 0;
      system.redistance();
      timer.registerTime("redistance");
      // }
      // lastRedistance++;
    } while (!system.advanceTime());

    if (advectionType == 0) {
      if (i % 5 == 0) {
        system.reseedParticles();
        timer.registerTime("reseed");

        system.adjustParticleRadius();
        timer.registerTime("radiusAdjust");
        std::cerr << " Final count: " << system.getParticleCount()
                  << " particles\n";
      }
      // system.printParticles();
    }

    system.printLevelSet();
    timer.registerTime("print");

    timer.evaluateComponentsTime();
    if (i % 5 == 0)
      timer.evaluateComponentsAverageTime();

    timer.logToFile("./results/log_weno.csv");
    Ramuh::FileWriter::writeArrayToFile("./results/particleCount.txt",
                                        particleHistory);
  }

  timer.evaluateComponentsAverageTime();

  return 0;
}
