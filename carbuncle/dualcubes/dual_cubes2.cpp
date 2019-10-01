#include <geometry/dual_marching.h>
#include <iostream>
#include "dual_cubes.h"

namespace Carbuncle {

// TODO: DualCubes2 massive tests!

DualCubes2::DualCubes2() {}

DualCubes2::DualCubes2(int size) {
  _h = 1. / size;
  _resolution = size;

  _cells.resize(size);
  _vertices.resize(size + 1);
  _ufaceNormals.resize(size + 1);
  _ufaceLocation.resize(size + 1);
  _vfaceNormals.resize(size);
  _vfaceLocation.resize(size);

  for (int i = 0; i < size; i++) {
    _cells[i].resize(size, 1e8);
    _vertices[i].resize(size + 1);
    _ufaceNormals[i].resize(size);
    _ufaceLocation[i].resize(size);
    _vfaceNormals[i].resize(size + 1);
    _vfaceLocation[i].resize(size + 1);
  }
  _vertices[size - 1].resize(size + 1);
  _ufaceNormals[size - 1].resize(size);
  _ufaceLocation[size - 1].resize(size);
}

void DualCubes2::initialize(Eigen::Array2d center, double radius) {
  // ====== Initialize cell surfaces
  for (int j = 0; j < _resolution; j++)
    for (int i = 0; i < _resolution; i++) {
      Eigen::Array2d h(_h, _h);
      Eigen::Array2d position = Eigen::Array2d(i, j).cwiseProduct(h) + h / 2.0;
      double distance = (position - center).matrix().norm() - radius;
      _cells[i][j] = std::min(_cells[i][j], distance);
    }

  // ====== Initialize edges and normals
  // Look for surface
  for (int j = 0; j < _resolution; j++)
    for (int i = 0; i < _resolution; i++) {
      int centerSign = (_cells[i][j] < 0) ? -1 : 1;
      if (i > 0) {
        int neighSign = (_cells[i - 1][j] < 0) ? -1 : 1;
        if (centerSign != neighSign) {
          // Compute intersection location
          double theta = _cells[i][j] - _cells[i - 1][j];
          double x = -_h * _cells[i - 1][j] / (theta) + (i - 0.5) * _h;
          double y = (j + 0.5) * _h;
          double denom = std::sqrt(x * x + y * y);
          _ufaceLocation[i][j] = Eigen::Array2d(x, y);
          _ufaceNormals[i][j] =
              Eigen::Vector2d(2 * (x - center[0]), 2 * (y - center[1]));
        }
      }
      if (i < _resolution - 1) {
        int neighSign = (_cells[i + 1][j] < 0) ? -1 : 1;
        if (centerSign != neighSign) {
          // Compute intersection location
          double theta = _cells[i + 1][j] - _cells[i][j];
          double x = -_h * _cells[i + 1][j] / (theta) + (i + 0.5) * _h;
          double y = (j + 0.5) * _h;
          double denom = std::sqrt(x * x + y * y);
          _ufaceLocation[i + 1][j] = Eigen::Array2d(x, y);
          _ufaceNormals[i + 1][j] =
              Eigen::Vector2d(2 * (x - center[0]), 2 * (y - center[1]));
        }
      }
      if (j > 0) {
        int neighSign = (_cells[i][j - 1] < 0) ? -1 : 1;
        if (centerSign != neighSign) {
          // Compute intersection location
          double theta = _cells[i][j] - _cells[i][j - 1];
          double x = (i + 0.5) * _h;
          double y = -_h * _cells[i][j - 1] / (theta) + (j - 0.5) * _h;
          double denom = std::sqrt(x * x + y * y);
          _vfaceLocation[i][j] = Eigen::Array2d(x, y);
          _vfaceNormals[i][j] =
              Eigen::Vector2d(2 * (x - center[0]), 2 * (y - center[1]));
        }
      }
      if (j < _resolution - 1) {
        int neighSign = (_cells[i][j + 1] < 0) ? -1 : 1;
        if (centerSign != neighSign) {
          // Compute intersection location
          double theta = _cells[i][j + 1] - _cells[i][j];
          double x = (i + 0.5) * _h;
          double y = -_h * _cells[i][j + 1] / (theta) + (j + 0.5) * _h;
          double denom = std::sqrt(x * x + y * y);
          _vfaceLocation[i][j + 1] = Eigen::Array2d(x, y);
          _vfaceNormals[i][j + 1] =
              Eigen::Vector2d(2 * (x - center[0]), 2 * (y - center[1]));
        }
      }
    }
}

void DualCubes2::printCells() {
  for (int j = 0; j < _resolution; j++) {
    for (int i = 0; i < _resolution; i++) {
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
  Ramuh::DualMarching2 surface(_resolution);
  std::cout << "======\n";

  for (int j = 0; j < _resolution - 1; j++) {
    for (int i = 0; i < _resolution - 1; i++) {
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
        auto x = surface.evaluateSquare(std::make_pair(i, j), normalPosition,
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