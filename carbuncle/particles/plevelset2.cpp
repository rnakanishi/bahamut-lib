#include <fluids/particle_levelset2.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

class PLevelSet : public Leviathan::ParticleLevelSet2 {
public:
  PLevelSet() : ParticleLevelSet2() {}

  PLevelSet(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain)
      : ParticleLevelSet2(gridSize, domain) {}

  void initializeLevelSet(Eigen::Array2d center, double radius) {
    auto &function = getCellScalarData(_phiId);
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = cellPosition(i);

      p -= center;
      function[i] = sqrt(p[0] * p[0] + p[1] * p[1]) - radius;
    }
  }

  void initializeGridVelocity() {
    auto &u = getFaceScalarData(0, _velocityId);
    auto &v = getFaceScalarData(1, _velocityId);

    // u velocities
    for (size_t i = 0; i < faceCount(0); i++) {
      auto p = facePosition(0, i);
      u[i] = p[1];

      // v velocities
      p = facePosition(1, i);
      v[i] = -p[0];
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

    for (size_t i = 0; i < particleCount(); i++) {
      auto pos = particlePosition(i);
      file << pos[0] << " " << pos[1] << " ";
      file << vel[i][0] << " " << vel[i][1] << " ";
      file << signal[i] << "\n";
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
      file << phi[i] << " ";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

protected:
};

int main(int argc, char const *argv[]) {
  PLevelSet system(Eigen::Array2i(64), Ramuh::BoundingBox2(-5, 5));

  system.initializeLevelSet(Eigen::Array2d(0, 1), 1.2);
  system.redistance();
  system.initializeGridVelocity();
  system.trackSurface();
  system.seedParticlesNearSurface();
  for (size_t i = 0; i < 10; i++) {
    system.interpolateVelocityToParticles();
    system.advectWeno();
    system.advectParticles();
    system.printParticles();
    system.printLevelSet();
  }

  return 0;
}
