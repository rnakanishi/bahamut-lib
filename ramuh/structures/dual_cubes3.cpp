#include <structures/dual_cubes.h>
#include <geometry/dual_marching.h>
#include <iostream>
#include <utility>

namespace Ramuh {

DualCubes3::DualCubes3() {}

DualCubes3::DualCubes3(int size) {
  _h = 1. / size;
  _resolution = size;

  _cells.resize(size);
  _vertices.resize(size + 1);
  _ufaceNormals.resize(size + 1);
  _ufaceLocation.resize(size + 1);
  _vfaceNormals.resize(size);
  _vfaceLocation.resize(size);
  _wfaceNormals.resize(size);
  _wfaceLocation.resize(size);

  for (int i = 0; i < size; i++) {
    _cells[i].resize(size);
    _vertices[i].resize(size + 1);
    _ufaceNormals[i].resize(size);
    _ufaceLocation[i].resize(size);
    _vfaceNormals[i].resize(size + 1);
    _vfaceLocation[i].resize(size + 1);
    _wfaceNormals[i].resize(size);
    _wfaceLocation[i].resize(size);

    _vertices[size - 1].resize(size + 1);
    _ufaceNormals[size - 1].resize(size);
    _ufaceLocation[size - 1].resize(size);
    for (int j = 0; j < size; j++) {
      _cells[i][j].resize(size, 1e8);
      _vertices[i][j].resize(size + 1);
      _ufaceNormals[i][j].resize(size);
      _ufaceLocation[i][j].resize(size);
      _vfaceNormals[i][j].resize(size);
      _vfaceLocation[i][j].resize(size);
      _wfaceNormals[i][j].resize(size + 1);
      _wfaceLocation[i][j].resize(size + 1);

      _vertices[i][size - 1].resize(size + 1);
      _ufaceNormals[size - 1][j].resize(size);
      _ufaceLocation[size - 1][j].resize(size);

      _vfaceNormals[i][size - 1].resize(size);
      _vfaceLocation[i][size - 1].resize(size);
      _vertices[size - 1][size - 1].resize(size + 1);
    }
  }
}

void DualCubes3::initialize(Eigen::Array3d center, double radius) {
  // ====== Initialize cell surfaces
  for (int k = 0; k < _resolution; k++)
    for (int j = 0; j < _resolution; j++)
      for (int i = 0; i < _resolution; i++) {
        Eigen::Array3d h(_h);
        Eigen::Array3d position =
            Eigen::Array3d(i, j, k).cwiseProduct(h) + h / 2.0;
        // double distance = (position - center).matrix().norm() - radius;

        // TORUS
        double x, y, z, x2, y2, z2;
        double distance;
        x = position[0];
        y = position[1];
        z = position[2];
        x2 = x * x;
        y2 = y * y;
        z2 = z * z;
        distance = ((x2 + y2) * (x2 + y2) - x2 + y2);
        distance = distance * distance + z2 - (1. / 100);

        _cells[i][j][k] = std::min(_cells[i][j][k], distance);
      }

  // ====== Initialize edges and normals
  // Look for surface
  for (int k = 0; k < _resolution; k++)
    for (int j = 0; j < _resolution; j++)
      for (int i = 0; i < _resolution; i++) {
        int centerSign = (_cells[i][j][k] < 0) ? -1 : 1;
        if (i > 0) {
          int neighSign = (_cells[i - 1][j][k] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j][k] - _cells[i - 1][j][k];
            double x = -_h * _cells[i - 1][j][k] / (theta) + (i - 0.5) * _h;
            double y = (j + 0.5) * _h;
            double z = (k + 0.5) * _h;
            _ufaceLocation[i][j][k] = Eigen::Array3d(x, y, z);

            // SPHERE
            // _ufaceNormals[i][j][k] = Eigen::Vector3d(
            // 2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));

            // TORUS
            _ufaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (4 * x * (x * x + y * y) - 2 * x) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * (4 * y * (x * x + y * y) + 2 * y) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * z);
          }
        }
        if (i < _resolution - 1) {
          int neighSign = (_cells[i + 1][j][k] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i + 1][j][k] - _cells[i][j][k];
            double x = -_h * _cells[i + 1][j][k] / (theta) + (i + 0.5) * _h;
            double y = (j + 0.5) * _h;
            double z = (k + 0.5) * _h;
            _ufaceLocation[i + 1][j][k] = Eigen::Array3d(x, y, z);

            // SPHERE
            // _ufaceNormals[i + 1][j][k] = Eigen::Vector3d(
            // 2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            _ufaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (4 * x * (x * x + y * y) - 2 * x) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * (4 * y * (x * x + y * y) + 2 * y) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * z);
          }
        }
        if (j > 0) {
          int neighSign = (_cells[i][j - 1][k] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j][k] - _cells[i][j - 1][k];
            double x = (i + 0.5) * _h;
            double y = -_h * _cells[i][j - 1][k] / (theta) + (j - 0.5) * _h;
            double z = (k + 0.5) * _h;
            _vfaceLocation[i][j][k] = Eigen::Array3d(x, y, z);
            // SPHERE
            // _vfaceNormals[i][j][k] = Eigen::Vector3d(
            // 2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            _vfaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (4 * x * (x * x + y * y) - 2 * x) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * (4 * y * (x * x + y * y) + 2 * y) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * z);
          }
        }
        if (j < _resolution - 1) {
          int neighSign = (_cells[i][j + 1][k] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j + 1][k] - _cells[i][j][k];
            double x = (i + 0.5) * _h;
            double y = -_h * _cells[i][j + 1][k] / (theta) + (j + 0.5) * _h;
            double z = (k + 0.5) * _h;
            _vfaceLocation[i][j + 1][k] = Eigen::Array3d(x, y, z);
            // SPHERE
            // _vfaceNormals[i][j + 1][k] = Eigen::Vector3d(
            // 2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            _vfaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (4 * x * (x * x + y * y) - 2 * x) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * (4 * y * (x * x + y * y) + 2 * y) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * z);
          }
        }
        if (k > 0) {
          int neighSign = (_cells[i][j][k - 1] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j][k] - _cells[i][j][k - 1];
            double x = (i + 0.5) * _h;
            double y = (j + 0.5) * _h;
            double z = -_h * _cells[i][j][k - 1] / (theta) + (k - 0.5) * _h;
            _wfaceLocation[i][j][k] = Eigen::Array3d(x, y, z);
            // SPHERE
            // _wfaceNormals[i][j][k] = Eigen::Vector3d(
            // 2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            _wfaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (4 * x * (x * x + y * y) - 2 * x) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * (4 * y * (x * x + y * y) + 2 * y) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * z);
          }
        }
        if (k < _resolution - 1) {
          int neighSign = (_cells[i][j][k + 1] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j][k + 1] - _cells[i][j][k];
            double x = (i + 0.5) * _h;
            double y = (j + 0.5) * _h;
            double z = -_h * _cells[i][j][k + 1] / (theta) + (k + 0.5) * _h;
            _wfaceLocation[i][j][k + 1] = Eigen::Array3d(x, y, z);
            // SPHERE
            // _wfaceNormals[i][j][k + 1] = Eigen::Vector3d(
            // 2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            _wfaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (4 * x * (x * x + y * y) - 2 * x) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * (4 * y * (x * x + y * y) + 2 * y) *
                    ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
                2 * z);
          }
        }
      }
}

void DualCubes3::printCells() {
  for (int k = 0; k < _resolution; k++) {
    for (int j = 0; j < _resolution; j++) {
      for (int i = 0; i < _resolution; i++) {
        if (!_ufaceNormals[i][j][k].isApprox(Eigen::Vector3d(0, 0, 0), 1e-8))
          std::cout << _ufaceLocation[i][j][k][0] << " "
                    << _ufaceLocation[i][j][k][1] << " "
                    << _ufaceLocation[i][j][k][2] << " "
                    << _ufaceNormals[i][j][k][0] << " "
                    << _ufaceNormals[i][j][k][1] << " "
                    << _ufaceNormals[i][j][k][2] << std::endl;
        if (!_vfaceNormals[i][j][k].isApprox(Eigen::Vector3d(0, 0, 0), 1e-8))
          std::cout << _vfaceLocation[i][j][k][0] << " "
                    << _vfaceLocation[i][j][k][1] << " "
                    << _vfaceLocation[i][j][k][2] << " "
                    << _vfaceNormals[i][j][k][0] << " "
                    << _vfaceNormals[i][j][k][1] << " "
                    << _vfaceNormals[i][j][k][2] << std::endl;
        if (!_wfaceNormals[i][j][k].isApprox(Eigen::Vector3d(0, 0, 0), 1e-8))
          std::cout << _wfaceLocation[i][j][k][0] << " "
                    << _wfaceLocation[i][j][k][1] << " "
                    << _wfaceLocation[i][j][k][2] << " "
                    << _wfaceNormals[i][j][k][0] << " "
                    << _wfaceNormals[i][j][k][1] << " "
                    << _wfaceNormals[i][j][k][2] << std::endl;
      }
    }
  }
}

void DualCubes3::extractSurface() {
  DualMarching3 surface(_resolution);
  std::cout << "======\n";

  for (int k = 0; k < _resolution - 1; k++) {
    for (int j = 0; j < _resolution - 1; j++) {
      for (int i = 0; i < _resolution - 1; i++) {
        std::vector<Eigen::Array3d> normalPosition;
        std::vector<Eigen::Vector3d> normal;
        bool isSurface = false;

        if (signChange(_cells[i][j][k], _cells[i + 1][j][k])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[i + 1][j][k]);
          normalPosition.emplace_back(_ufaceLocation[i + 1][j][k]);
        }
        if (signChange(_cells[i + 1][j + 1][k], _cells[i][j + 1][k])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[i + 1][j + 1][k]);
          normalPosition.emplace_back(_ufaceLocation[i + 1][j + 1][k]);
        }
        if (signChange(_cells[i][j][k + 1], _cells[i + 1][j][k + 1])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[i + 1][j][k + 1]);
          normalPosition.emplace_back(_ufaceLocation[i + 1][j][k + 1]);
        }
        if (signChange(_cells[i + 1][j + 1][k + 1], _cells[i][j + 1][k + 1])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[i + 1][j + 1][k + 1]);
          normalPosition.emplace_back(_ufaceLocation[i + 1][j + 1][k + 1]);
        }

        if (signChange(_cells[i + 1][j][k], _cells[i + 1][j + 1][k])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[i + 1][j + 1][k]);
          normalPosition.emplace_back(_vfaceLocation[i + 1][j + 1][k]);
        }
        if (signChange(_cells[i][j + 1][k], _cells[i][j][k])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[i][j + 1][k]);
          normalPosition.emplace_back(_vfaceLocation[i][j + 1][k]);
        }
        if (signChange(_cells[i + 1][j][k + 1], _cells[i + 1][j + 1][k + 1])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[i + 1][j + 1][k + 1]);
          normalPosition.emplace_back(_vfaceLocation[i + 1][j + 1][k + 1]);
        }
        if (signChange(_cells[i][j + 1][k + 1], _cells[i][j][k + 1])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[i][j + 1][k + 1]);
          normalPosition.emplace_back(_vfaceLocation[i][j + 1][k + 1]);
        }

        if (signChange(_cells[i][j][k], _cells[i][j][k + 1])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[i][j][k + 1]);
          normalPosition.emplace_back(_wfaceLocation[i][j][k + 1]);
        }
        if (signChange(_cells[i + 1][j][k], _cells[i + 1][j][k + 1])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[i + 1][j][k + 1]);
          normalPosition.emplace_back(_wfaceLocation[i + 1][j][k + 1]);
        }
        if (signChange(_cells[i + 1][j + 1][k], _cells[i + 1][j + 1][k + 1])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[i + 1][j + 1][k + 1]);
          normalPosition.emplace_back(_wfaceLocation[i + 1][j + 1][k + 1]);
        }
        if (signChange(_cells[i][j + 1][k], _cells[i][j + 1][k + 1])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[i][j + 1][k + 1]);
          normalPosition.emplace_back(_wfaceLocation[i][j + 1][k + 1]);
        }

        // Solve quadratic function
        if (isSurface) {
          auto x = surface.evaluateCube(std::make_tuple(i, j, k),
                                        normalPosition, normal);
          std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
        }
      }
    }
  }
  surface.reconstruct();
}

bool DualCubes3::signChange(double valueA, double valueB) {
  int signA = (valueA < 0) ? -1 : 1;
  int signB = (valueB < 0) ? -1 : 1;
  return signA != signB;
}
} // namespace Ramuh