#include <fluids/particle_levelset2.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utils/timer.hpp>

class PLevelSet : public Leviathan::ParticleLevelSet2 {
public:
  PLevelSet() : ParticleLevelSet2() {}

  PLevelSet(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain)
      : ParticleLevelSet2(gridSize, domain) {
    _dt = 1 / 30.;
  }

  void initializeLevelSet(Eigen::Array2d center, double radius) {
    auto &function = getCellScalarData(_phiId);
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getCellPosition(i);

      p -= center;
      // function[i] = -2;
      // if (sqrt(pow(p[0] - center[0], 2) + pow(p[1] - center[1], 2)) < radius)
      function[i] = sqrt(p[0] * p[0] + p[1] * p[1]) - radius;
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
    auto &u = getFaceScalarData(0, _velocityId);
    auto &v = getFaceScalarData(1, _velocityId);

    // u velocities
    for (size_t i = 0; i < faceCount(0); i++) {
      auto p = facePosition(0, i);
      u[i] = -p[1];

      // u[i] = 0;
    }

    // v velocities
    for (size_t j = 0; j < faceCount(1); j++) {
      auto p = facePosition(1, j);
      v[j] = p[0];
      // v[i] = -1;
    }
  }

  void printParticles() {
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/particles/2d/" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);
    auto &vel = getParticleArrayData(_particleVelocityId);
    auto &signal = getParticleScalarData(_particleSignalId);
    auto &radius = getParticleScalarData(_particleRadiusId);

    for (size_t i = 0; i < particleCount(); i++) {
      if (isActive(i)) {
        auto pos = getParticlePosition(i);
        file << pos[0] << " " << pos[1] << " ";
        file << vel[i][0] << " " << vel[i][1] << " ";
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
    filename << "results/particles/2d/ls" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);

    for (size_t i = 0; i < cellCount(); i++) {
      Eigen::Array2d pos = getCellPosition(i);
      file << pos[0] << " " << pos[1] << " " << phi[i] << "\n";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

  void defineVelocity(int i) {
    double time = _dt * i;
    auto &u = getFaceScalarData(0, _velocityId);
    auto &v = getFaceScalarData(1, _velocityId);

    // u velocities
    for (size_t i = 0; i < faceCount(0); i++) {
      auto p = facePosition(0, i);
      u[i] = -p[1];
      u[i] = pow(sin(M_PI * p[0]), 2) * sin(2 * M_PI * p[1]);
      u[i] *= cos(M_PI * time / 1.5);

      // u[i] = 0;
    }

    // v velocities
    for (size_t j = 0; j < faceCount(1); j++) {
      auto p = facePosition(1, j);
      v[j] = p[0];
      v[j] = -pow(sin(M_PI * p[1]), 2) * sin(2 * M_PI * p[0]);
      v[j] *= cos(M_PI * time / 1.5);
      // v[i] = -1;
    }
  }

protected:
};

int main(int argc, char const *argv[]) {
  // PLevelSet system(Eigen::Array2i(50), Ramuh::BoundingBox2(0, 1));
  PLevelSet system(Eigen::Array2i(50), Ramuh::BoundingBox2(-5, 5));

  // system.initializeLevelSet(Eigen::Array2d(.5, .75), 0.15);
  system.initializeLevelSet(Eigen::Array2d(0, 2), 1.5);
  system.redistance();
  system.computeCellsGradient();
  system.initializeGridVelocity();
  system.trackSurface();
  system.seedParticlesNearSurface();
  // system.attractParticles();
  system.printParticles();
  system.printLevelSet();

  Ramuh::Timer timer;

  int lastRedistance = 0;
  for (size_t i = 0; i <= 314; i++) {
    // system.defineVelocity(i);
    system.applyCfl();

    do {
      system.trackSurface();
      timer.registerTime("trackSurface");

      system.interpolateVelocityToParticles();
      timer.registerTime("interpolation");
      system.advectWeno();
      timer.registerTime("weno");
      system.advectParticles();
      timer.registerTime("particleAdvect");
    } while (!system.advanceTime());

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
      system.redistance();
    }
    lastRedistance++;

    system.printParticles();
    system.printLevelSet();
    timer.registerTime("print");

    system.reseedParticles();
    timer.registerTime("reseed");

    system.adjustParticleRadius();
    timer.registerTime("radiusAdjust");

    std::cerr << system.particleCount() << " particles\n";
    timer.evaluateComponentsTime();
  }

  return 0;
}
