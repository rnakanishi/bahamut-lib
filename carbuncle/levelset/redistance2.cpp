#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <fluids/levelset_fluid2.h>

class RedistanceClass : public Leviathan::LevelSetFluid2 {
public:
  RedistanceClass()
      : RedistanceClass(Eigen::Array2i(32, 32),
                        Ramuh::BoundingBox2::unitBox()) {}

  RedistanceClass(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain)
      : Leviathan::LevelSetFluid2(gridSize, domain) {}

  void initialize(Eigen::Array2d center, double radius) {
    auto &function = getScalarData(_phiId);
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getPosition(i);

      // p -= center;
      // function[i] = sqrt(p[0] * p[0] + p[1] * p[1]) - radius;

      // function[i] = p.abs().maxCoeff() - radius;
      function[i] = sqrt(p[0] * p[0] / (radius * radius) +
                         p[1] * p[1] / (radius * radius / 4)) -
                    1;
      function[i] *= 0.5 + (p[0] - center[0]) * (p[0] - center[0]) +
                     (p[1] - center[1]) * (p[1] - center[1]);
    }
  }

  void defineVelocity() {}

  void print() override {
    auto &phi = getScalarData(_phiId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/redistance/2d/" << count++;
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
  RedistanceClass cubes(
      Eigen::Array2i(20, 20),
      Ramuh::BoundingBox2(Eigen::Array2d(-5, -5), Eigen::Array2d(5, 5)));

  cubes.initialize(Eigen::Array2d(-1, 3), 4);
  // cubes.initialize(Eigen::Array2d(3.5, 2), 4);
  cubes.defineVelocity();

  cubes.computeCellsGradient();
  cubes.print();
  // return 1;

  // for (int i = 1; i <= 300; i++) {
  cubes.computeCellsGradient();
  cubes.redistance();

  cubes.print();
  // }

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
