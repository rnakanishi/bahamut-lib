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

    for (size_t i = 0; i < particleCount(); i++) {
      auto pos = particlePosition(i);
      file << pos[0] << " " << pos[1] << " ";
      file << vel[i][0] << " " << vel[i][1] << ";\n";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

protected:
};

int main(int argc, char const *argv[]) {
  PLevelSet system(Eigen::Array2i(32), Ramuh::BoundingBox2(-1, 1));

  system.seedParticles(Ramuh::BoundingBox2(-1, 1), 1000);
  system.initializeGridVelocity();
  for (size_t i = 0; i < 10; i++) {
    system.interpolateVelocityToParticles();
    system.advectParticles();
    system.printParticles();
  }

  return 0;
}
