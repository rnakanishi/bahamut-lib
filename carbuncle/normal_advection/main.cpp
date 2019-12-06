#include <iostream>
#include "level_set_reader.hpp"
#include "normal_particles.hpp"
#include <surface/dual_cubes.h>

int main(int argc, char const *argv[]) {
  /* code */

  Carbuncle::LevelSetReader reader(
      "/home/rnakanishi/Documents/meshes/cube3.sdf");
  Leviathan::LevelSetFluid3 levelset;

  reader.read(levelset);
  levelset.computeCellsGradient();

  std::cerr << "Seeding particles\n";
  Carbuncle::NormalParticles3 particles;
  particles.seedParticlesOverSurface(levelset);

  particles.print();

  // levelset.computeWenoGradient();
  // std::cerr << levelset.getDomain().getMin() <<
  // levelset.getDomain().getMax(); Leviathan::DualCubes cubes(levelset);
  // cubes.resetFileCounter();
  // cubes.setFolder("/home/rnakanishi/Documents/meshes/bunny/");
  // cubes.computeIntersectionAndNormals();
  // cubes.extractSurface();

  return 0;
}
