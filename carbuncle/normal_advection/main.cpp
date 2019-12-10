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

void printParticles(Carbuncle::NormalParticles3 particles) {
  std::ofstream file("matlab/particles.txt");
  auto &normalPair = particles.getPairMap();
  for (auto normal : normalPair) {
    auto origin = particles.getParticlePosition(normal.first);
    auto direction = particles.getParticlePosition(normal.second) - origin;
    file << origin.transpose() << " " << direction.transpose() << std::endl;
  }
  file.close();
}

void defineCellVelocity(Leviathan::LevelSetFluid3 levelset) {
  auto &u = levelset.getFaceArrayData(
      0, levelset.getFaceScalarLabelId("cellVelocity"));
  for (size_t faceId = 0; faceId < levelset.faceCount(0); faceId++) {
    auto ijk = levelset.idToijk(faceId);
    auto p = levelset.getFacePosition(0, faceId);
    u[faceId] = -p[1];
  }
  auto v = levelset.getFaceArrayData(
      0, levelset.getFaceScalarLabelId("cellVelocity"));
  for (size_t faceId = 0; faceId < levelset.faceCount(1); faceId++) {
    auto ijk = levelset.idToijk(faceId);
    auto p = levelset.getFacePosition(1, faceId);
    v[faceId] = p[0];
  }
}

void defineParticleVelocity(Carbuncle::NormalParticles3 particles) {
  auto vel = particles.getParticleArrayData("particleVelocity");
  for (size_t pid = 0; pid < particles.getParticleCount(); pid++) {
    if (particles.isActive(pid)) {
      auto p = particles.getParticlePosition(pid);
      vel[pid] = Eigen::Array3d(p[0], -p[1], 0);
    }
  }
}

int main(int argc, char const *argv[]) {
  Carbuncle::LevelSetReader reader(
      "/home/rnakanishi/Documents/meshes/zalesak_3d.sdf");
  Leviathan::LevelSetFluid3 levelset;

  reader.read(levelset);
  levelset.computeCellsGradient();

  // TODO: Change domain: insert dhe televelset to a predefined bounding box

  std::cerr << "Seeding particles\n";
  Carbuncle::NormalParticles3 particles;
  particles.seedParticlesOverSurface(levelset);

  particles.estimateCellNormals(levelset);
  Leviathan::DualCubes cubes(levelset);
  cubes.resetFileCounter();
  cubes.setFolder("/home/rnakanishi/Documents/meshes/bunny/");
  cubes.computeIntersectionAndNormals();
  cubes.extractSurface();

  printParticles(particles);
  printGradients(levelset);

  defineCellVelocity(levelset);
  defineParticleVelocity(particles);

  levelset.advectWeno();

  // TODO: Rotate levelset and particles
  // estimate surface normals with particles
  // Extract surface

  // TODO: Surface extraction should take particle position in consideration
  // when finding the zeroth levelset

  // levelset.computeWenoGradient();
  // std::cerr << levelset.getDomain().getMin() <<
  // levelset.getDomain().getMax();

  return 0;
}
