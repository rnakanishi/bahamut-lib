#include <geometry/dual_marching.h>
#include <iostream>
#include "dual_cubes.h"

namespace Carbuncle {

DualCubes3::DualCubes3()
    : DualCubes3(Eigen::Array3i(32, 32, 32), Ramuh::BoundingBox3::unitBox()) {}

DualCubes3::DualCubes3(Eigen::Array3i resolution)
    : DualCubes3(resolution, Ramuh::BoundingBox3::unitBox()) {}

DualCubes3::DualCubes3(Eigen::Array3i resolution, Ramuh::BoundingBox3 domain)
    : Leviathan::LevelSetFluid3(resolution, domain) {
  _h = _domain.size().cwiseQuotient(resolution.cast<double>());
  _resolution = resolution;
  newFaceArrayLabel("faceNormal");
  newFaceArrayLabel("facePosition");
}

void DualCubes3::initialize(Eigen::Array3d center, double radius) {
  this->initialize(center, radius, DualCubes3::ParametricSurface::SPHERE);
}

void DualCubes3::initialize(Eigen::Array3d center, double radius,
                            DualCubes3::ParametricSurface surface) {
  Eigen::Array3d domainMin = _domain.min();
  auto &_phi = getScalarVector("phi");

  // ====== Initialize cell surfaces
  for (int k = 0; k < _resolution[2]; k++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int i = 0; i < _resolution[0]; i++) {
        Eigen::Array3d h(_h);
        Eigen::Array3d position =
            domainMin + Eigen::Array3d(i, j, k).cwiseProduct(h) + h / 2.0;
        // double distance = (position - center).matrix().norm() - radius;

        double x, y, z, x2, y2, z2;
        double distance;
        x = position[0] - center[0];
        y = position[1] - center[1];
        z = position[2] - center[2];
        x2 = x * x;
        y2 = y * y;
        z2 = z * z;

        // TODO: Fix this initialization method. Extends this class into an
        // application and create another method to initialize properly
        switch (surface) {
        // SPHERE
        case DualCubes3::ParametricSurface::SPHERE:
          distance = x2 + y2 + z2 - radius * radius;
          break;
        case DualCubes3::ParametricSurface::TORUS:
          // TORUS
          distance = (1.5 * radius - std::sqrt(x2 + z2)) *
                         (1.5 * radius - std::sqrt(x2 + z2)) +
                     y2 - 0.75 * radius * radius;
          break;
        case DualCubes3::ParametricSurface::DOUBLE_TORUS:
          // DOUBLE TORUS
          distance = ((x2 + y2) * (x2 + y2) - x2 + y2);
          distance = distance * distance + z2 - (1. / 50);
          break;
        case DualCubes3::ParametricSurface::CUBE:
          // CUBE
          distance =
              std::max(std::max(std::fabs(x), std::fabs(y)), std::fabs(z)) -
              radius;

          break;
        default:
          distance = 1e8;
        }
        _phi[ijkToid(i, j, k)] = std::min(_phi[ijkToid(i, j, k)], distance);
      }
}

void DualCubes3::computeNormals() {
  auto &_phi = getScalarVector("phi");
  auto &_ufaceNormals = getFaceArrayVector(0, "faceNormal");
  auto &_vfaceNormals = getFaceArrayVector(1, "faceNormal");
  auto &_wfaceNormals = getFaceArrayVector(2, "faceNormal");
  for (int k = 0; k < _resolution[2]; k++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int i = 0; i < _resolution[0]; i++) {
        int centerSign = (_phi[ijkToid(i, j, k)] < 0) ? -1 : 1;
        double left, right, up, bot, back, front;
        if (i > 0) {
          int neighSign = (_phi[ijkToid(i - 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[0] =
                (_phi[ijkToid(i, j, k)] - _phi[ijkToid(i - 1, j, k)]) / _h[0];
            // Average first and then compute derivative for y
            // and z directions Y-direction
            up =
                (j < _resolution[1] - 1)
                    ? (_phi[ijkToid(i - 1, j + 1, k)] +
                       _phi[ijkToid(i, j + 1, k)]) /
                          2
                    : (_phi[ijkToid(i - 1, j, k)] + _phi[ijkToid(i, j, k)]) / 2;
            bot =
                (j > 0)
                    ? (_phi[ijkToid(i - 1, j - 1, k)] +
                       _phi[ijkToid(i, j - 1, k)]) /
                          2
                    : (_phi[ijkToid(i - 1, j, k)] + _phi[ijkToid(i, j, k)]) / 2;
            // Z-direction
            front =
                (k < _resolution[1] - 1)
                    ? (_phi[ijkToid(i - 1, j, k + 1)] +
                       _phi[ijkToid(i, j, k + 1)]) /
                          2
                    : (_phi[ijkToid(i - 1, j, k)] + _phi[ijkToid(i, j, k)]) / 2;
            back =
                (k > 0)
                    ? (_phi[ijkToid(i - 1, j, k - 1)] +
                       _phi[ijkToid(i, j, k - 1)]) /
                          2
                    : (_phi[ijkToid(i - 1, j, k)] + _phi[ijkToid(i, j, k)]) / 2;
            normal[1] = (up - bot) / (2 * _h[1]);
            normal[2] = (front - back) / (2 * _h[2]);
            _ufaceNormals[ijkToid(i, j, k)] = normal.normalized();
          }
        }
        if (i < _resolution[0] - 1) {
          int neighSign = (_phi[ijkToid(i + 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[0] =
                (_phi[ijkToid(i + 1, j, k)] - _phi[ijkToid(i, j, k)]) / _h[0];
            // Average first and then compute derivative for y
            // and z directions Y-direction
            up =
                (j < _resolution[1] - 1)
                    ? (_phi[ijkToid(i, j + 1, k)] +
                       _phi[ijkToid(i + 1, j + 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i + 1, j, k)]) / 2;
            bot =
                (j > 0)
                    ? (_phi[ijkToid(i, j - 1, k)] +
                       _phi[ijkToid(i + 1, j - 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i + 1, j, k)]) / 2;
            // Z-direction
            front =
                (k < _resolution[1] - 1)
                    ? (_phi[ijkToid(i, j, k + 1)] +
                       _phi[ijkToid(i + 1, j, k + 1)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i + 1, j, k)]) / 2;
            back =
                (k > 0)
                    ? (_phi[ijkToid(i, j, k - 1)] +
                       _phi[ijkToid(i + 1, j, k - 1)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i + 1, j, k)]) / 2;
            normal[1] = (up - bot) / (2 * _h[1]);
            normal[2] = (front - back) / (2 * _h[2]);
            _ufaceNormals[ijkToid(i, j, k)] = normal.normalized();
          }
        }
        if (j > 0) {
          int neighSign = (_phi[ijkToid(i, j - 1, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[1] =
                (_phi[ijkToid(i, j, k)] - _phi[ijkToid(i, j - 1, k)]) / _h[1];
            // Average first and then compute derivative for y
            // and z directions X-direction
            right =
                (i < _resolution[0] - 1)
                    ? (_phi[ijkToid(i + 1, j - 1, k)] +
                       _phi[ijkToid(i + 1, j, k)]) /
                          2
                    : (_phi[ijkToid(i, j - 1, k)] + _phi[ijkToid(i, j, k)]) / 2;
            left =
                (i > 0)
                    ? (_phi[ijkToid(i - 1, j - 1, k)] +
                       _phi[ijkToid(i - 1, j, k)]) /
                          2
                    : (_phi[ijkToid(i, j - 1, k)] + _phi[ijkToid(i, j, k)]) / 2;
            // Z-direction
            front =
                (k < _resolution[1] - 1)
                    ? (_phi[ijkToid(i, j - 1, k + 1)] +
                       _phi[ijkToid(i, j, k + 1)]) /
                          2
                    : (_phi[ijkToid(i, j - 1, k)] + _phi[ijkToid(i, j, k)]) / 2;
            back =
                (k > 0)
                    ? (_phi[ijkToid(i, j - 1, k - 1)] +
                       _phi[ijkToid(i, j, k - 1)]) /
                          2
                    : (_phi[ijkToid(i, j - 1, k)] + _phi[ijkToid(i, j, k)]) / 2;
            normal[0] = (right - left) / (2 * _h[0]);
            normal[2] = (front - back) / (2 * _h[2]);
            _vfaceNormals[ijkToid(i, j, k)] = normal.normalized();
          }
        }
        if (j < _resolution[1] - 1) {
          int neighSign = (_phi[ijkToid(i, j + 1, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[1] =
                (_phi[ijkToid(i, j + 1, k)] - _phi[ijkToid(i, j, k)]) / _h[1];
            // Average first and then compute derivative for y
            // and z directions Y-direction
            right =
                (i < _resolution[0] - 1)
                    ? (_phi[ijkToid(i + 1, j, k)] +
                       _phi[ijkToid(i + 1, j + 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i, j + 1, k)]) / 2;
            left =
                (i > 0)
                    ? (_phi[ijkToid(i - 1, j, k)] +
                       _phi[ijkToid(i - 1, j + 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i, j + 1, k)]) / 2;
            // Z-direction
            front =
                (k < _resolution[2] - 1)
                    ? (_phi[ijkToid(i, j, k + 1)] +
                       _phi[ijkToid(i, j + 1, k + 1)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i, j + 1, k)]) / 2;
            back =
                (k > 0)
                    ? (_phi[ijkToid(i, j, k - 1)] +
                       _phi[ijkToid(i, j + 1, k - 1)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i, j + 1, k)]) / 2;
            normal[0] = (right - left) / (2 * _h[0]);
            normal[2] = (front - back) / (2 * _h[2]);
            _vfaceNormals[ijkToid(i, j, k)] = normal.normalized();
          }
        }
        if (k > 0) {
          int neighSign = (_phi[ijkToid(i, j, k - 1)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[2] =
                (_phi[ijkToid(i, j, k)] - _phi[ijkToid(i, j, k - 1)]) / _h[2];
            // Average first and then compute derivative for y
            // and z directions X-direction
            right =
                (i < _resolution[0] - 1)
                    ? (_phi[ijkToid(i + 1, j - 1, k)] +
                       _phi[ijkToid(i + 1, j, k)]) /
                          2
                    : (_phi[ijkToid(i, j - 1, k)] + _phi[ijkToid(i, j, k)]) / 2;
            left =
                (i > 0)
                    ? (_phi[ijkToid(i - 1, j - 1, k)] +
                       _phi[ijkToid(i - 1, j, k)]) /
                          2
                    : (_phi[ijkToid(i, j - 1, k)] + _phi[ijkToid(i, j, k)]) / 2;
            // Y-direction
            up =
                (j < _resolution[1] - 1)
                    ? (_phi[ijkToid(i, j + 1, k - 1)] +
                       _phi[ijkToid(i, j + 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k - 1)] + _phi[ijkToid(i, j, k)]) / 2;
            bot =
                (j > 0)
                    ? (_phi[ijkToid(i, j - 1, k - 1)] +
                       _phi[ijkToid(i, j - 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k - 1)] + _phi[ijkToid(i, j, k)]) / 2;
            normal[0] = (right - left) / (2 * _h[0]);
            normal[1] = (up - bot) / (2 * _h[1]);
            _wfaceNormals[ijkToid(i, j, k)] = normal.normalized();
          }
        }
        if (k < _resolution[2] - 1) {
          int neighSign = (_phi[ijkToid(i, j, k + 1)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[1] =
                (_phi[ijkToid(i, j, k + 1)] - _phi[ijkToid(i, j, k)]) / _h[1];
            // Average first and then compute derivative for y
            // and z directions Y-direction
            right =
                (i < _resolution[0] - 1)
                    ? (_phi[ijkToid(i + 1, j, k)] +
                       _phi[ijkToid(i + 1, j + 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i, j + 1, k)]) / 2;
            left =
                (i > 0)
                    ? (_phi[ijkToid(i - 1, j, k)] +
                       _phi[ijkToid(i - 1, j + 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i, j + 1, k)]) / 2;
            // Y-direction
            up =
                (j < _resolution[1] - 1)
                    ? (_phi[ijkToid(i, j + 1, k + 1)] +
                       _phi[ijkToid(i, j + 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k + 1)] + _phi[ijkToid(i, j, k)]) / 2;
            bot =
                (j > 0)
                    ? (_phi[ijkToid(i, j - 1, k + 1)] +
                       _phi[ijkToid(i, j - 1, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k + 1)] + _phi[ijkToid(i, j, k)]) / 2;
            normal[0] = (right - left) / (2 * _h[0]);
            normal[1] = (up - bot) / (2 * _h[1]);
            _wfaceNormals[ijkToid(i, j, k)] = normal.normalized();
          }
        }
      }
}

void DualCubes3::computeIntersection() {
  auto &_phi = getScalarVector("phi");
  auto &_ufaceLocation = getFaceArrayVector(0, "facePosition");
  auto &_vfaceLocation = getFaceArrayVector(1, "facePosition");
  auto &_wfaceLocation = getFaceArrayVector(2, "facePosition");

  Eigen::Array3d domainMin = _domain.min();
  for (int k = 0; k < _resolution[2]; k++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int i = 0; i < _resolution[0]; i++) {
        int centerSign = (_phi[ijkToid(i, j, k)] < 0) ? -1 : 1;
        if (i > 0) {
          int neighSign = (_phi[ijkToid(i - 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i, j, k)] - _phi[ijkToid(i - 1, j, k)];
            double x = domainMin[0] +
                       -_h[0] * _phi[ijkToid(i - 1, j, k)] / (theta) +
                       (i - 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _ufaceLocation[ijkToid(i, j, k)] = Eigen::Array3d(x, y, z);
          }
        }
        if (i < _resolution[0] - 1) {
          int neighSign = (_phi[ijkToid(i + 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i + 1, j, k)] - _phi[ijkToid(i, j, k)];
            double x = domainMin[0] -
                       _h[0] * _phi[ijkToid(i + 1, j, k)] / (theta) +
                       (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _ufaceLocation[ijkToid(i + 1, j, k)] = Eigen::Array3d(x, y, z);
          }
        }
        if (j > 0) {
          int neighSign = (_phi[ijkToid(i, j - 1, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i, j, k)] - _phi[ijkToid(i, j - 1, k)];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] -
                       _h[1] * _phi[ijkToid(i, j - 1, k)] / (theta) +
                       (j - 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _vfaceLocation[ijkToid(i, j, k)] = Eigen::Array3d(x, y, z);
          }
        }
        if (j < _resolution[1] - 1) {
          int neighSign = (_phi[ijkToid(i, j + 1, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i, j + 1, k)] - _phi[ijkToid(i, j, k)];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] -
                       _h[1] * _phi[ijkToid(i, j + 1, k)] / (theta) +
                       (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _vfaceLocation[ijkToid(i, j + 1, k)] = Eigen::Array3d(x, y, z);
          }
        }
        if (k > 0) {
          int neighSign = (_phi[ijkToid(i, j, k - 1)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i, j, k)] - _phi[ijkToid(i, j, k - 1)];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] -
                       _h[2] * _phi[ijkToid(i, j, k - 1)] / (theta) +
                       (k - 0.5) * _h[2];
            _wfaceLocation[ijkToid(i, j, k)] = Eigen::Array3d(x, y, z);
          }
        }
        if (k < _resolution[2] - 1) {
          int neighSign = (_phi[ijkToid(i, j, k + 1)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i, j, k + 1)] - _phi[ijkToid(i, j, k)];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] +
                       -_h[2] * _phi[ijkToid(i, j, k + 1)] / (theta) +
                       (k + 0.5) * _h[2];
            _wfaceLocation[ijkToid(i, j, k + 1)] = Eigen::Array3d(x, y, z);
          }
        }
      }
}

void DualCubes3::printCells() {
  auto &_ufaceLocation = getFaceArrayVector(0, "facePosition");
  auto &_vfaceLocation = getFaceArrayVector(1, "facePosition");
  auto &_wfaceLocation = getFaceArrayVector(2, "facePosition");
  auto &_ufaceNormals = getFaceArrayVector(0, "faceNormal");
  auto &_vfaceNormals = getFaceArrayVector(1, "faceNormal");
  auto &_wfaceNormals = getFaceArrayVector(2, "faceNormal");

  for (int k = 0; k < _resolution[2]; k++) {
    for (int j = 0; j < _resolution[1]; j++) {
      for (int i = 0; i < _resolution[0]; i++) {
        if (!_ufaceNormals[ijkToid(i, j, k)].matrix().isApprox(
                Eigen::Vector3d(0, 0, 0), 1e-8))
          std::cout << _ufaceLocation[ijkToid(i, j, k)][0] << " "
                    << _ufaceLocation[ijkToid(i, j, k)][1] << " "
                    << _ufaceLocation[ijkToid(i, j, k)][2] << " "
                    << _ufaceNormals[ijkToid(i, j, k)][0] << " "
                    << _ufaceNormals[ijkToid(i, j, k)][1] << " "
                    << _ufaceNormals[ijkToid(i, j, k)][2] << std::endl;
        if (!_vfaceNormals[ijkToid(i, j, k)].matrix().isApprox(
                Eigen::Vector3d(0, 0, 0), 1e-8))
          std::cout << _vfaceLocation[ijkToid(i, j, k)][0] << " "
                    << _vfaceLocation[ijkToid(i, j, k)][1] << " "
                    << _vfaceLocation[ijkToid(i, j, k)][2] << " "
                    << _vfaceNormals[ijkToid(i, j, k)][0] << " "
                    << _vfaceNormals[ijkToid(i, j, k)][1] << " "
                    << _vfaceNormals[ijkToid(i, j, k)][2] << std::endl;
        if (!_wfaceNormals[ijkToid(i, j, k)].matrix().isApprox(
                Eigen::Vector3d(0, 0, 0), 1e-8))
          std::cout << _wfaceLocation[ijkToid(i, j, k)][0] << " "
                    << _wfaceLocation[ijkToid(i, j, k)][1] << " "
                    << _wfaceLocation[ijkToid(i, j, k)][2] << " "
                    << _wfaceNormals[ijkToid(i, j, k)][0] << " "
                    << _wfaceNormals[ijkToid(i, j, k)][1] << " "
                    << _wfaceNormals[ijkToid(i, j, k)][2] << std::endl;
      }
    }
  }
}

void DualCubes3::defineVelocity() {

  // u velocities
  for (int k = 0; k < _resolution[2]; k++) {
    for (int j = 0; j < _resolution[1]; j++) {
      for (int i = 0; i < _resolution[0] + 1; i++) {
        Eigen::Array3d facePosition;
        facePosition[1] = _domain.min()[1] + (j + 0.5) * _h[1];
        // _u[ijkToid(i, j, k)] = facePosition[1];
      }
    }
  }

  // v velocities
  for (int k = 0; k < _resolution[2]; k++) {
    for (int j = 0; j < _resolution[1] + 1; j++) {
      for (int i = 0; i < _resolution[0]; i++) {
        Eigen::Array3d facePosition;
        facePosition[0] = _domain.min()[0] + (i + 0.5) * _h[0];
        // _v[ijkToid(i, j, k)] = -facePosition[0];
      }
    }
  }

  // w velocities
  for (int k = 0; k < _resolution[2] + 1; k++) {
    for (int j = 0; j < _resolution[1]; j++) {
      for (int i = 0; i < _resolution[0]; i++) {
        Eigen::Array3d facePosition;
        // _w[ijkToid(i, j, k)] = 0;
      }
    }
  }
}

void DualCubes3::extractSurface() {
  Ramuh::DualMarching3 surface(_resolution);
  auto &_phi = getScalarVector("phi");
  auto &_ufaceLocation = getFaceArrayVector(0, "facePosition");
  auto &_vfaceLocation = getFaceArrayVector(1, "facePosition");
  auto &_wfaceLocation = getFaceArrayVector(2, "facePosition");
  auto &_ufaceNormals = getFaceArrayVector(0, "faceNormal");
  auto &_vfaceNormals = getFaceArrayVector(1, "faceNormal");
  auto &_wfaceNormals = getFaceArrayVector(2, "faceNormal");

  for (int k = 0; k < _resolution[2] - 1; k++) {
    for (int j = 0; j < _resolution[1] - 1; j++) {
      for (int i = 0; i < _resolution[0] - 1; i++) {
        std::vector<Eigen::Array3d> normalPosition;
        std::vector<Eigen::Vector3d> normal;
        bool isSurface = false;

        // yz-plane
        if (signChange(_phi[ijkToid(i, j, k)], _phi[ijkToid(i + 1, j, k)])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[ijkToid(i + 1, j, k)]);
          normalPosition.emplace_back(_ufaceLocation[ijkToid(i + 1, j, k)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j + 1, k)],
                       _phi[ijkToid(i, j + 1, k)])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[ijkToid(i + 1, j + 1, k)]);
          normalPosition.emplace_back(_ufaceLocation[ijkToid(i + 1, j + 1, k)]);
        }
        if (signChange(_phi[ijkToid(i, j, k + 1)],
                       _phi[ijkToid(i + 1, j, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[ijkToid(i + 1, j, k + 1)]);
          normalPosition.emplace_back(_ufaceLocation[ijkToid(i + 1, j, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j + 1, k + 1)],
                       _phi[ijkToid(i, j + 1, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[ijkToid(i + 1, j + 1, k + 1)]);
          normalPosition.emplace_back(
              _ufaceLocation[ijkToid(i + 1, j + 1, k + 1)]);
        }

        // xz-plane
        if (signChange(_phi[ijkToid(i + 1, j, k)],
                       _phi[ijkToid(i + 1, j + 1, k)])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[ijkToid(i + 1, j + 1, k)]);
          normalPosition.emplace_back(_vfaceLocation[ijkToid(i + 1, j + 1, k)]);
        }
        if (signChange(_phi[ijkToid(i, j + 1, k)], _phi[ijkToid(i, j, k)])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[ijkToid(i, j + 1, k)]);
          normalPosition.emplace_back(_vfaceLocation[ijkToid(i, j + 1, k)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j, k + 1)],
                       _phi[ijkToid(i + 1, j + 1, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[ijkToid(i + 1, j + 1, k + 1)]);
          normalPosition.emplace_back(
              _vfaceLocation[ijkToid(i + 1, j + 1, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i, j + 1, k + 1)],
                       _phi[ijkToid(i, j, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[ijkToid(i, j + 1, k + 1)]);
          normalPosition.emplace_back(_vfaceLocation[ijkToid(i, j + 1, k + 1)]);
        }

        // xy-plane
        if (signChange(_phi[ijkToid(i, j, k)], _phi[ijkToid(i, j, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[ijkToid(i, j, k + 1)]);
          normalPosition.emplace_back(_wfaceLocation[ijkToid(i, j, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j, k)],
                       _phi[ijkToid(i + 1, j, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[ijkToid(i + 1, j, k + 1)]);
          normalPosition.emplace_back(_wfaceLocation[ijkToid(i + 1, j, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j + 1, k)],
                       _phi[ijkToid(i + 1, j + 1, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[ijkToid(i + 1, j + 1, k + 1)]);
          normalPosition.emplace_back(
              _wfaceLocation[ijkToid(i + 1, j + 1, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i, j + 1, k)],
                       _phi[ijkToid(i, j + 1, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[ijkToid(i, j + 1, k + 1)]);
          normalPosition.emplace_back(_wfaceLocation[ijkToid(i, j + 1, k + 1)]);
        }
        // Solve quadratic function
        if (isSurface) {
          Eigen::Array3d cubeMin =
              _domain.min() +
              _h.cwiseProduct(Eigen::Array3d(i + 0.5, j + 0.5, k + 0.5));
          Eigen::Array3d cubeMax =
              _domain.min() +
              _h.cwiseProduct(Eigen::Array3d(i + 1.5, j + 1.5, k + 1.5));
          auto x = surface.evaluateCube(std::make_tuple(i, j, k),
                                        normalPosition, normal,
                                        Ramuh::BoundingBox3(cubeMin, cubeMax));
          // std::cout << x[0] << " " << x[1] << " " << x[2] <<
          // std::endl;
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
} // namespace Carbuncle