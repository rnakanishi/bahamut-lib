#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <structures/mac_grid1.h>
#include <fluids/levelset_fluid2.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

class Cip2 : public Leviathan::LevelSetFluid2 {
public:
  Cip2() : Cip2(Eigen::Array2i(32, 32), Ramuh::BoundingBox2::unitBox()) {}

  Cip2(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain)
      : Leviathan::LevelSetFluid2(gridSize, domain) {}

  void initialize(Eigen::Array2d center, double radius) {
    auto &phi = getScalarData(_phiId);
    auto &gradient = getArrayData(_gradientId);
    auto h = getH();
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getPosition(i);
      // phi[i] = 0.;
      // if (p[0] >= -4. && p[0] <= -3. && p[1] >= -4. && p[1] <= -3.)
      //   phi[i] = 2;

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

      phi[i] = std::sqrt(std::pow(p[0] - center[0], 2) +
                         std::pow(p[1] - center[1], 2)) -
               radius;
    }
    for (size_t id = 1; id < cellCount(); id++) {
      auto ij = idToij(id);
      int i = ij[0], j = ij[1];
      gradient[id][0] = (phi[id] - phi[ijToid(i - 1, j)]) / h[0];
      gradient[id][1] = (phi[id] - phi[ijToid(i, j - 1)]) / h[1];
    }

    defineVelocity();
  }

  void defineVelocity() {
    auto &uVelocity = getFaceScalarData(0, _velocityId);
    auto &vVelocity = getFaceScalarData(1, _velocityId);
    for (size_t i = 0; i < cellCount(); i++) {
      uVelocity[i] = 1.;
      vVelocity[i] = 1.;
    }
  }

  void print() {
    auto &phi = getScalarData(_phiId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/cip/2d/" << count++;
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
  Cip2 cubes(Eigen::Array2i(100), Ramuh::BoundingBox2(-5, 5));

  cubes.initialize(Eigen::Array2d(0, 0), 1);

  cubes.print();

  for (int i = 1; i <= 100; i++) {
    cubes.advectCip();

    cubes.print();
  }

  //   auto surface = cubes.marc hingTetrahedra();
  //   Ramuh::FileWriter writer;
  //   std::ostringstream objname;
  //   objname << "results/marching/tetra_" << std::setfill('0') << std::setw(4)
  //   << 0
  //   << ".obj";
  //   writer.writeMeshModel(surface, objname.str());
  // }
  return 0;
}
