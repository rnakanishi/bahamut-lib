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

  initializeCube(cubes, Eigen::Array2d(0, 0), 3, ParametricSurface::SQUARE);
  initializeGradientsAtIntersection(cubes, Eigen::Array2d(0, 0), 3,
                                    ParametricSurface::SQUARE);
  defineVelocity(cubes);
  cubes.setFolder("results/dualSquares/");
  cubes.computeCellsGradient();

  Carbuncle::NormalParticles2 particles;
  particles.seedParticlesOverSurface(cubes);

  printParticles(particles);

  // cubes.redistance();

  // cubes.computeCellsGradient();
  // cubes.computeIntersectionAndNormals();

  cubes.extractSurface();
  return 1;

  for (int i = 1; i <= 15; i++) {
    cubes.advectWeno();

    if (i % 10 == 0)
      cubes.redistance();
    cubes.computeWenoGradient();
    cubes.computeIntersectionAndNormals();
    cubes.extractSurface();
  }

  return 0;
}
