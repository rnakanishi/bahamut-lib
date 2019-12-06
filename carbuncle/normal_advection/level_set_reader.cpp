#include "level_set_reader.hpp"
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <geometry/bounding_box.h>

namespace Carbuncle {

LevelSetReader::LevelSetReader(std::string filename) { _filename = filename; }

void LevelSetReader::read(Leviathan::LevelSetFluid3 &levelset) {
  std::ifstream file(_filename);

  if (!file.is_open()) {
    std::cerr << "File not opened: " << _filename << std::endl;
    return;
  }
  std::string line;
  Eigen::Array3i dimension;
  Eigen::Array3d origin;
  double dx;
  file >> dimension[0] >> dimension[1] >> dimension[2];
  file >> origin[0] >> origin[1] >> origin[2];
  file >> dx;

  // Read configuration
  std::cout << "Dimensions: " << dimension.matrix().transpose() << std::endl;
  std::cout << "Origin: " << origin.matrix().transpose() << std::endl;
  std::cout << "dx: " << dx << std::endl;

  // Initialize levelset with given configuration
  levelset = Leviathan::LevelSetFluid3(
      dimension,
      Ramuh::BoundingBox3(origin, origin + dx * dimension.cast<double>()));

  // Read data
  for (size_t k = 0; k < dimension[2]; k++) {
    for (size_t j = 0; j < dimension[1]; j++) {
      for (size_t i = 0; i < dimension[0]; i++) {
        double value;
        file >> value;
        levelset.setPhiValue(Eigen::Array3i(i, j, k), value);
      }
    }
  }
  file.close();
  std::cerr << "SDF file " << _filename << " read completely\n";
}

} // namespace Carbuncle