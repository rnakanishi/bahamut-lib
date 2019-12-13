#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <iomanip>
#include "dual_cubes.h"
#include <surface/dual_squares.h>
#include <structures/line_mesh.hpp>
#include "initialization.hpp"
#include "tangent_particles2.hpp"

void printParticles(Carbuncle::TangentParticles2 particles) {
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
      Eigen::Array2i(40, 40),
      Ramuh::BoundingBox2(Eigen::Array2d(-5, -5), Eigen::Array2d(5, 5)));
  Carbuncle::TangentParticles2 particles(cubes);

  Eigen::Array2d center = Eigen::Array2d(0, 0);
  double radius = 1.5;

  initializeCube(cubes, center, radius, ParametricSurface::SQUARE);
  initializeGradientsAtIntersection(cubes, center, radius,
                                    ParametricSurface::SQUARE);
  cubes.setFolder("results/dualSquares/");

  defineCellsVelocity(cubes);

  particles.seedParticlesOverSurface(cubes);
  defineParticlesVelocity(particles);
  printParticles(particles);

  // cubes.computeCellsGradient();
  // cubes.computeIntersectionAndNormals();

  cubes.extractSurface();
  particles.extractSurface(cubes);

  for (int i = 1; i <= 1; i++) {
    cubes.advectWeno();
    particles.advectParticles();

    cubes.computeIntersectionAndNormals();
    particles.fixLevelsetGradients(cubes);

    // After levelset and particles advection, cell gradient should be corrected
    // using particle information
    cubes.computeIntersectionAndNormals();
    // Correct computed gradients with particle normals
    cubes.computeCellsGradient();
    cubes.computeIntersectionAndNormals();
    // particles.estimateCellNormals(cubes);
    cubes.extractSurface();
    particles.extractSurface(cubes);
  }
  cubes.computeCellsGradient();
  cubes.computeIntersectionAndNormals();
  particles.fixLevelsetGradients(cubes);
  cubes.extractSurface();

  particles.extractSurface(cubes);
  return 0;
}
