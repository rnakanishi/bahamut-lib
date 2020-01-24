#include <fluids/particle_levelset3.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

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

  void setMeshFolder(std::string folder) { meshFolder = folder; }

  std::string getMeshFolder() { return meshFolder; }

protected:
  int lscount;
  int pfileCount;
  std::string baseFolder, meshFolder;
};