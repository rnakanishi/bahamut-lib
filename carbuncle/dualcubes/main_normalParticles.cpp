#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <iomanip>
#include "dual_cubes.h"
#include <surface/dual_squares.h>
#include "initialization.hpp"
#include "normal_particles2.hpp"

void printParticles(Carbuncle::NormalParticles2 particles) {
  std::ofstream file("matlab/particles.txt");
  auto &normalPair = particles.getPairMap();
  for (auto normal : normalPair) {
    auto origin = particles.getParticlePosition(normal.first);
    auto direction = particles.getParticlePosition(normal.second) - origin;
    file << origin.transpose() << " " << direction.transpose() << std::endl;
  }
  file.close();
}

int main(int argc, char const *argv[]) {
  Leviathan::DualSquares cubes(
      Eigen::Array2i(25, 25),
      Ramuh::BoundingBox2(Eigen::Array2d(-5, -5), Eigen::Array2d(5, 5)));
  Carbuncle::NormalParticles2 particles;

  initializeCube(cubes, Eigen::Array2d(0, 0), 3, ParametricSurface::SQUARE);
  initializeGradientsAtIntersection(cubes, Eigen::Array2d(0, 0), 3,
                                    ParametricSurface::SQUARE);
  cubes.setFolder("results/dualSquares/");

  defineCellsVelocity(cubes);

  particles.seedParticlesOverSurface(cubes);
  defineParticlesVelocity(particles);
  printParticles(particles);

  // cubes.computeCellsGradient();
  // cubes.computeIntersectionAndNormals();

  cubes.extractSurface();

  for (int i = 1; i <= 50; i++) {
    cubes.advectWeno();
    particles.advectParticles();

    // After levelset and particles advection, cell gradient should be corrected
    // using particle information
    cubes.computeIntersectionAndNormals();
    // Correct computed gradients with particle normals
  }
  printParticles(particles);
  cubes.computeCellsGradient();
  cubes.computeIntersectionAndNormals();
  cubes.extractSurface();

  return 0;
}
