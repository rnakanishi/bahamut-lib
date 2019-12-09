#include <iostream>
#include "level_set_reader.hpp"
#include "normal_particles.hpp"
#include <surface/dual_cubes.h>
#include <fstream>

void printGradients(Leviathan::LevelSetFluid3 &levelset) {
  auto gradients = levelset.getCellArrayData("cellGradient");
  auto h = levelset.getH();
  auto domain = levelset.getDomain();
  std::ofstream file("matlab/gradients");

  for (size_t cellId = 0; cellId < levelset.cellCount(); cellId++) {
    auto position = levelset.getCellPosition(cellId);
    file << position.transpose() << " " << gradients[cellId].transpose();
    file << std::endl;
  }

  file.close();
}

int main(int argc, char const *argv[]) {
  /* code */

  Carbuncle::LevelSetReader reader(
      "/home/rnakanishi/Documents/meshes/icosphere.sdf");
  Leviathan::LevelSetFluid3 levelset;

  reader.read(levelset);
  levelset.computeCellsGradient();

  std::cerr << "Seeding particles\n";
  Carbuncle::NormalParticles3 particles;
  particles.seedParticlesOverSurface(levelset);

  particles.print();
  printGradients(levelset);

  // levelset.computeWenoGradient();
  // std::cerr << levelset.getDomain().getMin() <<
  // levelset.getDomain().getMax(); Leviathan::DualCubes cubes(levelset);
  // cubes.resetFileCounter();
  // cubes.setFolder("/home/rnakanishi/Documents/meshes/bunny/");
  // cubes.computeIntersectionAndNormals();
  // cubes.extractSurface();

  return 0;
}
