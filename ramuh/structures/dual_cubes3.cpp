#include <structures/dual_cubes.h>
#include <geometry/dual_marching.h>
#include <iostream>

namespace Ramuh {

DualCubes3::DualCubes3()
    : DualCubes3(Eigen::Array3i(32, 32, 32), BoundingBox3::unitBox()) {}

DualCubes3::DualCubes3(Eigen::Array3i resolution)
    : DualCubes3(resolution, BoundingBox3::unitBox()) {}

DualCubes3::DualCubes3(Eigen::Array3i resolution, BoundingBox3 domain) {
  _domain = domain;
  _h = _domain.size().cwiseQuotient(resolution.cast<double>());
  _resolution = resolution;

  _cells.resize(resolution[0]);
  _vertices.resize(resolution[0] + 1);
  _ufaceNormals.resize(resolution[0] + 1);
  _ufaceLocation.resize(resolution[0] + 1);
  _vfaceNormals.resize(resolution[0]);
  _vfaceLocation.resize(resolution[0]);
  _wfaceNormals.resize(resolution[0]);
  _wfaceLocation.resize(resolution[0]);

  for (int i = 0; i < resolution[0]; i++) {
    _cells[i].resize(resolution[1]);
    _vertices[i].resize(resolution[1] + 1);
    _ufaceNormals[i].resize(resolution[1]);
    _ufaceLocation[i].resize(resolution[1]);
    _vfaceNormals[i].resize(resolution[1] + 1);
    _vfaceLocation[i].resize(resolution[1] + 1);
    _wfaceNormals[i].resize(resolution[1]);
    _wfaceLocation[i].resize(resolution[1]);

    _vertices[resolution[0] - 1].resize(resolution[1] + 1);
    _ufaceNormals[resolution[0] - 1].resize(resolution[1]);
    _ufaceLocation[resolution[0] - 1].resize(resolution[1]);
    for (int j = 0; j < resolution[1]; j++) {
      _cells[i][j].resize(resolution[2], 1e8);
      _vertices[i][j].resize(resolution[2] + 1);
      _ufaceNormals[i][j].resize(resolution[2]);
      _ufaceLocation[i][j].resize(resolution[2]);
      _vfaceNormals[i][j].resize(resolution[2]);
      _vfaceLocation[i][j].resize(resolution[2]);
      _wfaceNormals[i][j].resize(resolution[2] + 1);
      _wfaceLocation[i][j].resize(resolution[2] + 1);

      _vertices[i][resolution[1] - 1].resize(resolution[2] + 1);
      _ufaceNormals[resolution[1] - 1][j].resize(resolution[2]);
      _ufaceLocation[resolution[1] - 1][j].resize(resolution[2]);

      _vfaceNormals[i][resolution[1] - 1].resize(resolution[2]);
      _vfaceLocation[i][resolution[1] - 1].resize(resolution[2]);
      _vertices[resolution[1] - 1][resolution[1] - 1].resize(resolution[2] + 1);
    }
  }
}

void DualCubes3::initialize(Eigen::Array3d center, double radius) {
  Eigen::Array3d domainMin = _domain.min();
  // ====== Initialize cell surfaces
  for (int k = 0; k < _resolution[2]; k++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int i = 0; i < _resolution[0]; i++) {
        Eigen::Array3d h(_h);
        Eigen::Array3d position =
            domainMin + Eigen::Array3d(i, j, k).cwiseProduct(h) + h / 2.0;
        // double distance = (position - center).matrix().norm() - radius;

        // TORUS
        double x, y, z, x2, y2, z2;
        double distance;
        x = position[0] - center[0];
        y = position[1] - center[1];
        z = position[2] - center[2];
        x2 = x * x;
        y2 = y * y;
        z2 = z * z;

        // SPHERE
        distance = x2 + y2 + z2 - radius * radius;

        // TORUS
        // distance = ((x2 + y2) * (x2 + y2) - x2 + y2);
        // distance = distance * distance + z2 - (1. / 100);

        _cells[i][j][k] = std::min(_cells[i][j][k], distance);
      }

  // ====== Initialize edges and normals
  // Look for surface
  for (int k = 0; k < _resolution[2]; k++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int i = 0; i < _resolution[0]; i++) {
        int centerSign = (_cells[i][j][k] < 0) ? -1 : 1;
        if (i > 0) {
          int neighSign = (_cells[i - 1][j][k] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j][k] - _cells[i - 1][j][k];
            double x = domainMin[0] + -_h[0] * _cells[i - 1][j][k] / (theta) +
                       (i - 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _ufaceLocation[i][j][k] = Eigen::Array3d(x, y, z);

            // SPHERE
            _ufaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));

            // TORUS
            // _ufaceNormals[i][j][k] = Eigen::Vector3d(
            //     2 * (4 * x * (x * x + y * y) - 2 * x) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * (4 * y * (x * x + y * y) + 2 * y) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * z);
          }
        }
        if (i < _resolution[0] - 1) {
          int neighSign = (_cells[i + 1][j][k] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i + 1][j][k] - _cells[i][j][k];
            double x = domainMin[0] - _h[0] * _cells[i + 1][j][k] / (theta) +
                       (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _ufaceLocation[i + 1][j][k] = Eigen::Array3d(x, y, z);

            // SPHERE
            _ufaceNormals[i + 1][j][k] = Eigen::Vector3d(
                2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            // _ufaceNormals[i][j][k] = Eigen::Vector3d(
            //     2 * (4 * x * (x * x + y * y) - 2 * x) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * (4 * y * (x * x + y * y) + 2 * y) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * z);
          }
        }
        if (j > 0) {
          int neighSign = (_cells[i][j - 1][k] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j][k] - _cells[i][j - 1][k];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] - _h[1] * _cells[i][j - 1][k] / (theta) +
                       (j - 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _vfaceLocation[i][j][k] = Eigen::Array3d(x, y, z);
            // SPHERE
            _vfaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            // _vfaceNormals[i][j][k] = Eigen::Vector3d(
            //     2 * (4 * x * (x * x + y * y) - 2 * x) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * (4 * y * (x * x + y * y) + 2 * y) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * z);
          }
        }
        if (j < _resolution[1] - 1) {
          int neighSign = (_cells[i][j + 1][k] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j + 1][k] - _cells[i][j][k];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] - _h[1] * _cells[i][j + 1][k] / (theta) +
                       (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _vfaceLocation[i][j + 1][k] = Eigen::Array3d(x, y, z);
            // SPHERE
            _vfaceNormals[i][j + 1][k] = Eigen::Vector3d(
                2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            // _vfaceNormals[i][j][k] = Eigen::Vector3d(
            //     2 * (4 * x * (x * x + y * y) - 2 * x) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * (4 * y * (x * x + y * y) + 2 * y) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * z);
          }
        }
        if (k > 0) {
          int neighSign = (_cells[i][j][k - 1] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j][k] - _cells[i][j][k - 1];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] - _h[2] * _cells[i][j][k - 1] / (theta) +
                       (k - 0.5) * _h[2];
            _wfaceLocation[i][j][k] = Eigen::Array3d(x, y, z);
            // SPHERE
            _wfaceNormals[i][j][k] = Eigen::Vector3d(
                2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            // _wfaceNormals[i][j][k] = Eigen::Vector3d(
            //     2 * (4 * x * (x * x + y * y) - 2 * x) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * (4 * y * (x * x + y * y) + 2 * y) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * z);
          }
        }
        if (k < _resolution[2] - 1) {
          int neighSign = (_cells[i][j][k + 1] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _cells[i][j][k + 1] - _cells[i][j][k];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + -_h[2] * _cells[i][j][k + 1] / (theta) +
                       (k + 0.5) * _h[2];
            _wfaceLocation[i][j][k + 1] = Eigen::Array3d(x, y, z);
            // SPHERE
            _wfaceNormals[i][j][k + 1] = Eigen::Vector3d(
                2 * (x - center[0]), 2 * (y - center[1]), 2 * (z - center[2]));
            // TORUS
            // _wfaceNormals[i][j][k] = Eigen::Vector3d(
            //     2 * (4 * x * (x * x + y * y) - 2 * x) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * (4 * y * (x * x + y * y) + 2 * y) *
            //         ((x * x + y * y) * (x * x + y * y) - x * x + y * y),
            //     2 * z);
          }
        }
      }
}

void DualCubes3::printCells() {
  for (int k = 0; k < _resolution[2]; k++) {
    for (int j = 0; j < _resolution[1]; j++) {
      for (int i = 0; i < _resolution[0]; i++) {
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

  for (int k = 0; k < _resolution[2] - 1; k++) {
    for (int j = 0; j < _resolution[1] - 1; j++) {
      for (int i = 0; i < _resolution[0] - 1; i++) {
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