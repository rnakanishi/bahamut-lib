#include <geometry/dual_marching.h>
#include <iostream>
#include "dual_cubes.h"

namespace Carbuncle {

// TODO: DualCubes2 massive tests!

DualCubes2::DualCubes2()
    : DualCubes2(Eigen::Array2i(32, 32), Ramuh::BoundingBox2::unitBox()) {}

DualCubes2::DualCubes2(Eigen::Array2i resolution)
    : DualCubes2(resolution, Ramuh::BoundingBox2::unitBox()) {}

DualCubes2::DualCubes2(Eigen::Array2i resolution, Ramuh::BoundingBox2 domain)
    : Leviathan::LevelSetFluid2(resolution, domain) {
  _h = _domain.getSize().cwiseQuotient(resolution.cast<double>());
  newFaceArrayLabel("faceNormal");
  newFaceArrayLabel("facePosition");
}

void DualCubes2::initialize(Eigen::Array2d center, double radius,
                            DualCubes2::ParametricSurface surface) {
  Eigen::Array2d domainMin = _domain.getMin();
  auto &_phi = getCellScalarData("phi");

  // ====== Initialize cell surfaces
  for (int j = 0; j < _gridSize[1]; j++)
    for (int i = 0; i < _gridSize[0]; i++) {
      Eigen::Array2d h(_h);
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
      case DualCubes2::ParametricSurface::CIRCLE:
        distance = x2 + y2 - radius * radius;
        break;
      case DualCubes2::ParametricSurface::SQUARE:
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
      _phi[ijToid(i, j)] = std::min(_phi[ijToid(i, j)], distance);
    }
}

void DualCubes2::computeIntersectionAndNormals() {
  auto &phi = getCellScalarData("phi");
  auto &gradient = getCellArrayData("cellGradient");
  auto &ufaceLocation = getFaceArrayData(0, "facePosition");
  auto &vfaceLocation = getFaceArrayData(1, "facePosition");
  auto &ufaceNormals = getFaceArrayData(0, "faceNormal");
  auto &vfaceNormals = getFaceArrayData(1, "faceNormal");

  Eigen::Array2d domainMin = _domain.getMin();
  // #pragma omp parallel for
  for (int j = 0; j < _gridSize[1]; j++)
    for (int i = 0; i < _gridSize[0]; i++) {
      int centerSign = (phi[ijToid(i, j)] < 0) ? -1 : 1;
      int id = ijToid(i, j);
      if (i < _gridSize[0] - 1) {
        int neighSign = (phi[ijToid(i + 1, j)] < 0) ? -1 : 1;
        if (centerSign != neighSign) {
          // Compute intersection location
          double theta = phi[ijToid(i + 1, j)] - phi[ijToid(i, j)];
          double x = domainMin[0] - _h[0] * phi[ijToid(i, j)] / (theta) +
                     (i + 0.5) * _h[0];
          double y = domainMin[1] + (j + 0.5) * _h[1];
          ufaceLocation[faceijToid(0, i + 1, j)] = Eigen::Array2d(x, y);
          ufaceNormals[faceijToid(0, i + 1, j)] =
              (gradient[id] - gradient[ijToid(i + 1, j)]) * phi[id] / theta +
              gradient[id];
        }
      }
      if (j < _gridSize[1] - 1) {
        int neighSign = (phi[ijToid(i, j + 1)] < 0) ? -1 : 1;
        if (centerSign != neighSign) {
          // Compute intersection location
          double theta = phi[ijToid(i, j + 1)] - phi[ijToid(i, j)];
          double x = domainMin[0] + (i + 0.5) * _h[0];
          double y = domainMin[1] - _h[1] * phi[ijToid(i, j)] / (theta) +
                     (j + 0.5) * _h[1];
          vfaceLocation[faceijToid(1, i, j + 1)] = Eigen::Array2d(x, y);
          vfaceNormals[faceijToid(1, i, j + 1)] =
              (gradient[id] - gradient[ijToid(i, j + 1)]) * phi[id] / theta +
              gradient[id];
        }
      }
    }
}

void DualCubes2::printCells() {
  for (int j = 0; j < _gridSize[1]; j++) {
    for (int i = 0; i < _gridSize[0]; i++) {
      if (!_ufaceNormals[i][j].isApprox(Eigen::Vector2d(0, 0), 1e-8))
        std::cout << _ufaceLocation[i][j][0] << " " << _ufaceLocation[i][j][1]
                  << " " << _ufaceNormals[i][j][0] << " "
                  << _ufaceNormals[i][j][1] << std::endl;
      if (!_vfaceNormals[i][j].isApprox(Eigen::Vector2d(0, 0), 1e-8))
        std::cout << _vfaceLocation[i][j][0] << " " << _vfaceLocation[i][j][1]
                  << " " << _vfaceNormals[i][j][0] << " "
                  << _vfaceNormals[i][j][1] << std::endl;
    }
  }
}

void DualCubes2::extractSurface() {
  Ramuh::DualMarching2 surface(_gridSize);
  std::cout << "======\n";

  for (int j = 0; j < _gridSize[1] - 1; j++) {
    for (int i = 0; i < _gridSize[0] - 1; i++) {
      std::vector<Eigen::Array2d> normalPosition;
      std::vector<Eigen::Vector2d> normal;
      bool isSurface = false;

      if (signChange(_cells[i][j], _cells[i + 1][j])) {
        isSurface = true;
        normal.emplace_back(_ufaceNormals[i + 1][j]);
        normalPosition.emplace_back(_ufaceLocation[i + 1][j]);
      }
      if (signChange(_cells[i + 1][j], _cells[i + 1][j + 1])) {
        isSurface = true;
        normal.emplace_back(_vfaceNormals[i + 1][j + 1]);
        normalPosition.emplace_back(_vfaceLocation[i + 1][j + 1]);
      }
      if (signChange(_cells[i + 1][j + 1], _cells[i][j + 1])) {
        isSurface = true;
        normal.emplace_back(_ufaceNormals[i + 1][j + 1]);
        normalPosition.emplace_back(_ufaceLocation[i + 1][j + 1]);
      }
      if (signChange(_cells[i][j + 1], _cells[i][j])) {
        isSurface = true;
        normal.emplace_back(_vfaceNormals[i][j + 1]);
        normalPosition.emplace_back(_vfaceLocation[i][j + 1]);
      }

      // Solve quadratic function
      if (isSurface) {
        auto x = surface.evaluateSquare(Eigen::Array2i(i, j), normalPosition,
                                        normal);
        std::cout << x[0] << " " << x[1] << std::endl;
      }
    }
  }
}

bool DualCubes2::signChange(double valueA, double valueB) {
  int signA = (valueA < 0) ? -1 : 1;
  int signB = (valueB < 0) ? -1 : 1;
  return signA != signB;
}
} // namespace Carbuncle