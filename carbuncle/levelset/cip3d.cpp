#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <structures/mac_grid1.h>
#include <fluids/levelset_fluid3.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

class Cip3 : public Leviathan::LevelSetFluid3 {
public:
  Cip3() : Cip3(Eigen::Array3i(32), Ramuh::BoundingBox3::unitBox()) {}

  Cip3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain)
      : Leviathan::LevelSetFluid3(gridSize, domain) {
    _originalDt = _dt = 2 * M_PI / 628;
  }

  void initialize(Eigen::Array3d center, double radius) {
    auto &phi = getScalarData(_phiId);
    auto &gradient = getArrayData(_gradientId);
    auto h = getH();
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getPosition(i);
      phi[i] = -2;
      // if (p[0] >= -4. && p[0] <= -3. && p[1] >= -4. && p[1] <= -3.)
      if (pow(p[0] - center[0], 2) + pow(p[1] - center[1], 2) +
              pow(p[2] - center[2], 2) <
          radius * radius)
        phi[i] = 5;
      double s = radius / 4.;
      double s6 = radius / 6.;
      if (p[0] >= center[0] - s && p[0] <= center[0] + s &&
          p[1] < center[1] + 2 * s6 && p[1] > center[1] - radius - s6 &&
          p[2] >= center[2] - radius && p[2] <= center[2] + radius)
        phi[i] = -2;

      // double A = 2;
      // double sigma[2] = {0.5, 0.5};
      // phi[i] =
      //     A * std::exp(
      //             -(std::pow(p[0] - center[0], 2) / (2 * sigma[0] * sigma[0])
      //             +
      //               std::pow(p[1] - center[1], 2) / (2 * sigma[1] *
      //               sigma[1])));
      // if (phi[i] < 1e-16)
      //   phi[i] = 0;

      // phi[i] = std::sqrt(std::pow(p[0] - center[0], 2) +
      //                    std::pow(p[1] - center[1], 2)) -
      //          radius;
    }
    for (size_t id = 1; id < cellCount(); id++) {
      auto ijk = idToijk(id);
      int i = ijk[0], j = ijk[1], k = ijk[2];
      gradient[id][0] = (phi[id] - phi[ijkToid(i - 1, j, k)]) / h[0];
      gradient[id][1] = (phi[id] - phi[ijkToid(i, j - 1, k)]) / h[1];
      gradient[id][2] = (phi[id] - phi[ijkToid(i, j, k - 1)]) / h[2];
    }

    defineVelocity();
  }

  void defineVelocity() {
    auto &uVelocity = getFaceScalarData(0, _velocityId);
    auto &vVelocity = getFaceScalarData(1, _velocityId);
    auto &wVelocity = getFaceScalarData(2, _velocityId);
    for (size_t id = 0; id < faceCount(0); id++) {
      auto p = facePosition(0, id);
      uVelocity[id] = -p[1];
    }

    for (size_t id = 0; id < faceCount(1); id++) {
      auto p = facePosition(1, id);
      vVelocity[id] = p[0];
    }

    for (size_t id = 0; id < faceCount(1); id++) {
      auto p = facePosition(2, id);
      wVelocity[id] = 0;
    }
  }

  void print() {
    auto &phi = getScalarData(_phiId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/cip/3d/" << count++;
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
  Cip3 cubes(Eigen::Array3i(80), Ramuh::BoundingBox3(-5, 5));

  cubes.initialize(Eigen::Array3d(0, 2, 0), 2);
  cubes.computeCellsGradient();
  cubes.print();
  // cubes.redistance();

  for (int i = 1; i <= 15; i++) {
    cubes.applyCfl();
    do {
      cubes.advectCip();

      // cubes.advectWeno();
      // if (i % 10 == 0)
      //   cubes.redistance();

    } while (!cubes.advanceTime());
    // if (i % 10 == 0)
    //   cubes.redistance();

    cubes.print();
  }

  return 0;
}