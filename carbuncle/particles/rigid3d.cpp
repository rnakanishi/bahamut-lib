#include <fluids/particle_levelset3.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utils/timer.hpp>

class PLevelSet3 : public Leviathan::ParticleLevelSet3 {
public:
  PLevelSet3() : ParticleLevelSet3() {}

  PLevelSet3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain)
      : ParticleLevelSet3(gridSize, domain) {
    // _dt = 1 / 30.;
  }

  void initializeLevelSet(Eigen::Array3d center, double radius) {
    auto &function = getCellScalarData(_phiId);
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getCellPosition(i);

      p -= center;
      // function[i] = -2;
      // if (sqrt(pow(p[0] - center[0], 2) + pow(p[1] - center[1], 2)) < radius)
      function[i] = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - radius;
      // function[i] = 5;

      double s = radius / 4.;
      if ((p[0] >= -s && p[0] <= s) && (p[1] <= s && p[1] >= -radius)) {
        function[i] =
            std::min(abs(s - p[1]), std::min(abs(p[0] + s), abs(s - p[0])));
      }
      if (function[i] < 0) {
        if (p[1] <= s && p[1] >= -radius)
          function[i] = -std::min(abs(function[i]),
                                  std::min(abs(p[0] + s), abs(s - p[0])));
        if (p[0] >= -s && p[0] <= s)
          function[i] = -std::min(abs(function[i]), abs(s - p[1]));
      }
    }
  }

  void initializeGridVelocity() {
    auto &u = getFaceScalarData(0, _cellVelocityId);
    auto &v = getFaceScalarData(1, _cellVelocityId);
    auto &w = getFaceScalarData(2, _velocityId);

    // u velocities
    for (size_t i = 0; i < faceCount(0); i++) {
      auto p = facePosition(0, i);
      u[i] = p[1];
    }

    // v velocities
    for (size_t j = 0; j < faceCount(1); j++) {
      auto p = facePosition(1, j);
      v[j] = -p[0];
    }

    for (size_t k = 0; k < faceCount(1); k++) {
      auto p = facePosition(1, j);
      w[k] = 0;
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

    for (size_t i = 0; i < getParticleCount(); i++) {
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
    double time = _dt * i;
    auto &u = getFaceScalarData(0, _cellVelocityId);
    auto &v = getFaceScalarData(1, _cellVelocityId);
    auto &w = getFaceScalarData(2, _cellVelocityId);

    // u velocities
    for (size_t i = 0; i < faceCount(0); i++) {
      auto p = getFacePosition(0, i);
      u[i] = p[1];
      u[i] = 2 * pow(sin(M_PI * p[0]), 2) * sin(2 * M_PI * p[1]) *
             sin(2 * M_PI * p[2]);
      u[i] *= cos(M_PI * time / 3.0);

      // u[i] = 0;
    }

    // v velocities
    for (size_t j = 0; j < faceCount(1); j++) {
      auto p = getFacePosition(1, j);
      v[j] = -p[0];
      v[j] = -pow(sin(M_PI * p[1]), 2) * sin(2 * M_PI * p[0]) *
             sin(2 * M_PI * p[2]);
      v[j] *= cos(M_PI * time / 3.0);
      // v[i] = -1;
    }

    // w velocities
    for (size_t k = 0; k < faceCount(1); k++) {
      auto p = getFacePosition(2, k);
      w[k] = -p[0];
      w[k] = -pow(sin(M_PI * p[2]), 2) * sin(2 * M_PI * p[0]) *
             sin(2 * M_PI * p[1]);
      w[k] *= cos(M_PI * time / 3.0);
      // v[i] = -1;
    }
  }
  void setBaseFolder(std::string folder) { baseFolder = folder; }

protected:
  std::string baseFolder;
};

int main(int argc, char const *argv[]) {
  PLevelSet3 system(Eigen::Array3i(150), Ramuh::BoundingBox3(-5, 5));
  system.setBaseFolder("results/particles/3d/pls_rigid");

  system.initializeLevelSet(Eigen::Array3d(0, 2), 1.5);
  system.redistance();
  system.computeCellsGradient();
  system.initializeGridVelocity();
  system.trackSurface();
  system.seedParticlesNearSurface();
  // system.attractParticles();
  system.printParticles();
  system.printLevelSet();
  system.printGradients();
  system.printVelocities();

  Ramuh::Timer timer;

  int lastRedistance = 0;
  for (size_t i = 0; i <= 500; i++) {
    // system.defineVelocity(i);
    system.applyCfl();

    do {
      system.trackSurface();
      timer.registerTime("trackSurface");

      system.interpolateVelocityToParticles();
      timer.registerTime("interpolation");
      // system.advectSemiLagrangian();
      system.advectWeno();
      timer.registerTime("cellAdvection");
      system.advectParticles();
      timer.registerTime("particleAdvect");

      if (system.correctLevelSetWithParticles()) {
        std::cerr << "......................................... CORRECTED\n";
        timer.registerTime("correction");

        lastRedistance = 0;
        system.redistance();
        timer.registerTime("redistance");

        system.correctLevelSetWithParticles();
        timer.registerTime("correction2");
      } else
        timer.registerTime("correction");

      if (lastRedistance >= 5) {
        lastRedistance = 0;
        std::cerr << "Redistance iteration\n";
        system.redistance();
      }
      lastRedistance++;

      system.reseedParticles();
      timer.registerTime("reseed");

      system.adjustParticleRadius();
      timer.registerTime("radiusAdjust");
    } while (!system.advanceTime());
    system.printParticles();
    system.printLevelSet();
    system.printGradients();
    system.printVelocities();
    timer.registerTime("print");

    std::cerr << system.getParticleCount() << " particles\n";
    timer.evaluateComponentsTime();
  }

  return 0;
}
