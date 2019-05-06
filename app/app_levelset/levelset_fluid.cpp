#include "levelset_fluid.h"
#include <fstream>

LevelSetFluid2::LevelSetFluid2() : Ramuh::LevelSet2() {}

LevelSetFluid3::LevelSetFluid3() : Ramuh::LevelSet3() {}

void LevelSetFluid3::writeVelocityField() {
  std::ofstream file;
  file.open("results/velocityField");
  if (file.is_open()) {
    for (int id = 0; id < cellCount(); id++) {
      Eigen::Array3i ijk = idToijk(id);
      int i, j, k;
      i = ijk[0];
      j = ijk[1];
      k = ijk[2];

      file << (_u[i][j][k] + _u[i + 1][j][k]) / 2 << " ";
      file << (_v[i][j][k] + _v[i][j + 1][k]) / 2 << " ";
      file << (_w[i][j][k] + _u[i][j][k + 1]) / 2 << "\n";
    }
  }
}