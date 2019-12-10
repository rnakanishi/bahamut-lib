#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <iomanip>
#include "dual_cubes.h"
#include <surface/dual_squares.h>

enum class ParametricSurface : int { CIRCLE, SQUARE };

void defineVelocity(Leviathan::DualSquares &cubes) {}

void initializeCube(Leviathan::DualSquares &squares, Eigen::Array2d center,
                    double radius, ParametricSurface surface) {

  Eigen::Array2d domainMin = squares.getDomain().getMin();
  auto &_phi = squares.getCellScalarData("phi");
  auto gridSize = squares.getGridSize();
  auto h = squares.getH();

  // ====== Initialize cell surfaces
  for (int j = 0; j < gridSize[1]; j++)
    for (int i = 0; i < gridSize[0]; i++) {
      Eigen::Array2d position =
          domainMin + Eigen::Array2d(i, j).cwiseProduct(h) + h / 2.0;
      // double distance = (position - center).matrix().norm() - radius;

      double x, y, x2, y2;
      double distance;
      x = position[0] - center[0];
      y = position[1] - center[1];
      x2 = x * x;
      y2 = y * y;

      double r2;
      double r;
      // TODO: Fix this initialization method. Extends this class into an
      // application and create another method to initialize properly
      switch (surface) {
      // SPHERE
      case ParametricSurface::CIRCLE:
        distance = x2 + y2 - radius * radius;
        break;
      case ParametricSurface::SQUARE:
        // CUBE
        distance = std::max(std::fabs(x), std::fabs(y)) - radius;
        if (distance > 0) {
          position = position.abs();
          distance = 0.0;
          x = std::max(0.0, position[0] - radius);
          y = std::max(0.0, position[1] - radius);
          distance = sqrt(x * x + y * y);
        }
        break;
      default:
        distance = 1e8;
      }
      _phi[squares.ijToid(i, j)] =
          std::min(_phi[squares.ijToid(i, j)], distance);
    }
}

int main(int argc, char const *argv[]) {
  Leviathan::DualSquares cubes(
      Eigen::Array2i(10, 10),
      Ramuh::BoundingBox2(Eigen::Array2d(-5, -5), Eigen::Array2d(5, 5)));

  initializeCube(cubes, Eigen::Array2d(0, 0), 3, ParametricSurface::SQUARE);

  defineVelocity(cubes);
  cubes.setFolder("results/dualSquares/");

  // cubes.redistance();

  cubes.computeCellsGradient();
  cubes.computeIntersectionAndNormals();

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
