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
      p -= center;

      // function[i] = p.abs().maxCoeff() - radius;
      function[i] = sqrt(p[0] * p[0] + p[1] * p[1]) - radius * radius;
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
      Eigen::Array2i(32, 32),
      Ramuh::BoundingBox2(Eigen::Array2d(-2, -2), Eigen::Array2d(2, 2)));

  cubes.initialize(Eigen::Array2d(0, 0), 1.0);
  cubes.defineVelocity();

  cubes.computeCellsGradient();
  cubes.print();
  // return 1;

  for (int i = 1; i <= 100; i++) {
    cubes.computeCellsGradient();
    cubes.redistance();

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
