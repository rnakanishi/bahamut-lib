#include <geometry/dual_marching.h>
#include <omp.h>
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
  newArrayLabel("cellGradient");
  newFaceArrayLabel("faceNormal");
  newFaceArrayLabel("facePosition");
}

void DualCubes3::initialize(Eigen::Array3d center, double radius) {
  this->initialize(center, radius, DualCubes3::ParametricSurface::SPHERE);
}

void DualCubes3::initialize(Eigen::Array3d center, double radius,
                            DualCubes3::ParametricSurface surface) {
  Eigen::Array3d domainMin = _domain.min();
  auto &_phi = getScalarData("phi");

  // ====== Initialize cell surfaces
  for (int k = 0; k < _gridSize[2]; k++)
    for (int j = 0; j < _gridSize[1]; j++)
      for (int i = 0; i < _gridSize[0]; i++) {
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
          if (distance > 0) {
            position = position.abs();
            distance = 0.0;
            x = std::max(0.0, position[0] - radius);
            y = std::max(0.0, position[1] - radius);
            z = std::max(0.0, position[2] - radius);
            distance = sqrt(x * x + y * y + z * z);
          }
          break;
        default:
          distance = 1e8;
        }
        _phi[ijkToid(i, j, k)] = std::min(_phi[ijkToid(i, j, k)], distance);
      }
}

void DualCubes3::analyticNormals(Eigen::Array3d center, double radius,
                                 DualCubes3::ParametricSurface surface) {
  Eigen::Array3d domainMin = _domain.min();
  auto &phi = getScalarData("phi");
  auto &ufacePosition = getFaceArrayData(0, "facePosition");
  auto &vfacePosition = getFaceArrayData(1, "facePosition");
  auto &wfacePosition = getFaceArrayData(2, "facePosition");
  auto &ufaceNormal = getFaceArrayData(0, "faceNormal");
  auto &vfaceNormal = getFaceArrayData(1, "faceNormal");
  auto &wfaceNormal = getFaceArrayData(2, "faceNormal");

  for (int k = 1; k < _gridSize[2] - 1; k++)
    for (int j = 1; j < _gridSize[1] - 1; j++)
      for (int i = 1; i < _gridSize[0] - 1; i++) {
        Eigen::Array3d h(_h);
        Eigen::Array3d position =
            domainMin + Eigen::Array3d(i, j, k).cwiseProduct(h) + h / 2.0;
        int centerSign = (phi[ijkToid(i, j, k)] < 0) ? -1 : 1;
        Eigen::Array3d normal;
        Eigen::Array3d intersec;
        int neighSign;
        switch (surface) {
        // SPHERE
        case DualCubes3::ParametricSurface::SPHERE:
          break;
        case DualCubes3::ParametricSurface::TORUS:
          break;
        case DualCubes3::ParametricSurface::DOUBLE_TORUS:
          break;
        case DualCubes3::ParametricSurface::CUBE:
          neighSign = (phi[ijkToid(i - 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            intersec = facePosition(0, i, j, k);
            double theta = phi[ijkToid(i, j, k)] - phi[ijkToid(i - 1, j, k)];
            intersec[0] = domainMin[0] +
                          -_h[0] * phi[ijkToid(i - 1, j, k)] / (theta) +
                          (i - 0.5) * _h[0];
            ufacePosition[faceijkToid(0, i, j, k)] = intersec;
            normal = Eigen::Array3d(0);
            if (ufacePosition[faceijkToid(0, i, j, k)][0] < center[0])
              normal[0] = -1;
            else
              normal[0] = 1;
            ufaceNormal[faceijkToid(0, i, j, k)] = normal;
          }

          neighSign = (phi[ijkToid(i, j - 1, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            intersec = facePosition(1, i, j, k);
            double theta = phi[ijkToid(i, j, k)] - phi[ijkToid(i, j - 1, k)];
            intersec[1] = domainMin[1] +
                          -_h[1] * phi[ijkToid(i, j - 1, k)] / (theta) +
                          (j - 0.5) * _h[1];
            vfacePosition[faceijkToid(1, i, j, k)] = intersec;
            normal = Eigen::Array3d(0);
            if (vfacePosition[faceijkToid(1, i, j, k)][1] < center[1])
              normal[1] = -1;
            else
              normal[1] = 1;
            vfaceNormal[faceijkToid(1, i, j, k)] = normal;
          }

          neighSign = (phi[ijkToid(i, j, k - 1)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            intersec = facePosition(2, i, j, k);
            double theta = phi[ijkToid(i, j, k)] - phi[ijkToid(i, j, k - 1)];
            intersec[2] = domainMin[2] +
                          -_h[2] * phi[ijkToid(i, j, k - 1)] / (theta) +
                          (k - 0.5) * _h[2];
            wfacePosition[faceijkToid(2, i, j, k)] = intersec;
            normal = Eigen::Array3d(0);
            if (wfacePosition[faceijkToid(2, i, j, k)][2] < center[2])
              normal[2] = -1;
            else
              normal[2] = 1;
            wfaceNormal[faceijkToid(2, i, j, k)] = normal;
          }
          break;
        default:
          break;
        }
      }
}

void DualCubes3::computeNormals() {
  auto &_phi = getScalarData("phi");
  auto &_ufaceNormals = getFaceArrayData(0, "faceNormal");
  auto &_vfaceNormals = getFaceArrayData(1, "faceNormal");
  auto &_wfaceNormals = getFaceArrayData(2, "faceNormal");
  // TODO: Store normal location only on negative levelset cells
  // TODO: Voundary treatment for discrete normal computation
  for (int k = 0; k < _gridSize[2]; k++)
    for (int j = 0; j < _gridSize[1]; j++)
      for (int i = 0; i < _gridSize[0]; i++) {
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
                (j < _gridSize[1] - 1)
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
                (k < _gridSize[1] - 1)
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
            _ufaceNormals[faceijkToid(0, i, j, k)] = normal.normalized();
          }
        }
        if (i < _gridSize[0] - 1) {
          int neighSign = (_phi[ijkToid(i + 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[0] =
                (_phi[ijkToid(i + 1, j, k)] - _phi[ijkToid(i, j, k)]) / _h[0];
            // Average first and then compute derivative for y
            // and z directions Y-direction
            up =
                (j < _gridSize[1] - 1)
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
                (k < _gridSize[1] - 1)
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
            _ufaceNormals[faceijkToid(0, i, j, k)] = normal.normalized();
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
                (i < _gridSize[0] - 1)
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
                (k < _gridSize[1] - 1)
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
            _vfaceNormals[faceijkToid(1, i, j, k)] = normal.normalized();
          }
        }
        if (j < _gridSize[1] - 1) {
          int neighSign = (_phi[ijkToid(i, j + 1, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[1] =
                (_phi[ijkToid(i, j + 1, k)] - _phi[ijkToid(i, j, k)]) / _h[1];
            // Average first and then compute derivative for y
            // and z directions Y-direction
            right =
                (i < _gridSize[0] - 1)
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
                (k < _gridSize[2] - 1)
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
            _vfaceNormals[faceijkToid(1, i, j, k)] = normal.normalized();
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
                (i < _gridSize[0] - 1)
                    ? (_phi[ijkToid(i + 1, j, k - 1)] +
                       _phi[ijkToid(i + 1, j, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k - 1)] + _phi[ijkToid(i, j, k)]) / 2;
            left =
                (i > 0)
                    ? (_phi[ijkToid(i - 1, j, k - 1)] +
                       _phi[ijkToid(i - 1, j, k)]) /
                          2
                    : (_phi[ijkToid(i, j, k - 1)] + _phi[ijkToid(i, j, k)]) / 2;
            // Y-direction
            up =
                (j < _gridSize[1] - 1)
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
            _wfaceNormals[faceijkToid(2, i, j, k)] = normal.normalized();
          }
        }
        if (k < _gridSize[2] - 1) {
          int neighSign = (_phi[ijkToid(i, j, k + 1)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            Eigen::Vector3d normal;
            normal[1] =
                (_phi[ijkToid(i, j, k + 1)] - _phi[ijkToid(i, j, k)]) / _h[1];
            // Average first and then compute derivative for y
            // and z directions Y-direction
            right =
                (i < _gridSize[0] - 1)
                    ? (_phi[ijkToid(i + 1, j, k)] +
                       _phi[ijkToid(i + 1, j, k + 1)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i, j, k + 1)]) / 2;
            left =
                (i > 0)
                    ? (_phi[ijkToid(i - 1, j, k)] +
                       _phi[ijkToid(i - 1, j, k + 1)]) /
                          2
                    : (_phi[ijkToid(i, j, k)] + _phi[ijkToid(i, j, k + 1)]) / 2;
            // Y-direction
            up =
                (j < _gridSize[1] - 1)
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
            _wfaceNormals[faceijkToid(2, i, j, k)] = normal.normalized();
          }
        }
      }
}

void DualCubes3::computeCellsGradient() {
  auto &phi = getScalarData("phi");
  auto &gradient = getArrayData("cellGradient");
  auto h = getH();
#pragma omp parallel for
  for (int k = 0; k < _gridSize[2]; k++)
    for (int j = 0; j < _gridSize[1]; j++)
      for (int i = 0; i < _gridSize[0]; i++) {
        int id = ijkToid(i, j, k);
        if (i < _gridSize[0])
          gradient[id][0] = (phi[ijkToid(i + 1, j, k)] - phi[id]) / h[0];
        else
          gradient[id] = (phi[id] - phi[ijkToid(i - 1, j, k)]) / h[0];

        if (j < _gridSize[1])
          gradient[id][1] = (phi[ijkToid(i, j + 1, k)] - phi[id]) / h[1];
        else
          gradient[id] = (phi[id] - phi[ijkToid(i, j - 1, k)]) / h[1];

        if (k < _gridSize[2])
          gradient[id][2] = (phi[ijkToid(i, j, k + 1)] - phi[id]) / h[2];
        else
          gradient[id] = (phi[id] - phi[ijkToid(i, j, k - 1)]) / h[2];
      }
}

void DualCubes3::computeIntersectionAndNormals() {
  auto &phi = getScalarData("phi");
  auto &gradient = getArrayData("cellGradient");
  auto &ufaceLocation = getFaceArrayData(0, "facePosition");
  auto &vfaceLocation = getFaceArrayData(1, "facePosition");
  auto &wfaceLocation = getFaceArrayData(2, "facePosition");
  auto &ufaceNormals = getFaceArrayData(0, "faceNormal");
  auto &vfaceNormals = getFaceArrayData(1, "faceNormal");
  auto &wfaceNormals = getFaceArrayData(2, "faceNormal");

  Eigen::Array3d domainMin = _domain.min();
#pragma omp parallel for
  for (int k = 0; k < _gridSize[2]; k++)
    for (int j = 0; j < _gridSize[1]; j++)
      for (int i = 0; i < _gridSize[0]; i++) {
        int centerSign = (phi[ijkToid(i, j, k)] < 0) ? -1 : 1;
        int id = ijkToid(i, j, k);
        if (i < _gridSize[0] - 1) {
          int neighSign = (phi[ijkToid(i + 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = phi[ijkToid(i + 1, j, k)] - phi[ijkToid(i, j, k)];
            double x = domainMin[0] - _h[0] * phi[ijkToid(i, j, k)] / (theta) +
                       (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            ufaceLocation[faceijkToid(0, i + 1, j, k)] =
                Eigen::Array3d(x, y, z);
            ufaceNormals[faceijkToid(0, i + 1, j, k)] =
                (gradient[id] - gradient[ijkToid(i + 1, j, k)]) * phi[id] /
                    theta +
                gradient[id];
          }
        }
        if (j < _gridSize[1] - 1) {
          int neighSign = (phi[ijkToid(i, j + 1, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = phi[ijkToid(i, j + 1, k)] - phi[ijkToid(i, j, k)];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] - _h[1] * phi[ijkToid(i, j, k)] / (theta) +
                       (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            vfaceLocation[faceijkToid(1, i, j + 1, k)] =
                Eigen::Array3d(x, y, z);
            vfaceNormals[faceijkToid(1, i, j + 1, k)] =
                (gradient[id] - gradient[ijkToid(i, j + 1, k)]) * phi[id] /
                    theta +
                gradient[id];
          }
        }
        if (k < _gridSize[2] - 1) {
          int neighSign = (phi[ijkToid(i, j, k + 1)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = phi[ijkToid(i, j, k + 1)] - phi[ijkToid(i, j, k)];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + -_h[2] * phi[ijkToid(i, j, k)] / (theta) +
                       (k + 0.5) * _h[2];
            wfaceLocation[faceijkToid(2, i, j, k + 1)] =
                Eigen::Array3d(x, y, z);
            wfaceNormals[faceijkToid(2, i, j, k + 1)] =
                (gradient[id] - gradient[ijkToid(i, j, k + 1)]) * phi[id] /
                    theta +
                gradient[id];
          }
        }
      }
}

void DualCubes3::computeIntersection() {
  auto &_phi = getScalarData("phi");
  auto &_ufaceLocation = getFaceArrayData(0, "facePosition");
  auto &_vfaceLocation = getFaceArrayData(1, "facePosition");
  auto &_wfaceLocation = getFaceArrayData(2, "facePosition");

  Eigen::Array3d domainMin = _domain.min();
  for (int k = 0; k < _gridSize[2]; k++)
    for (int j = 0; j < _gridSize[1]; j++)
      for (int i = 0; i < _gridSize[0]; i++) {
        int centerSign = (_phi[ijkToid(i, j, k)] < 0) ? -1 : 1;
        if (i > 0) {
          int neighSign = (_phi[ijkToid(i - 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i, j, k)] - _phi[ijkToid(i - 1, j, k)];
            double x = domainMin[0] -
                       _h[0] * _phi[ijkToid(i - 1, j, k)] / (theta) +
                       (i - 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _ufaceLocation[faceijkToid(0, i, j, k)] = Eigen::Array3d(x, y, z);
          }
        }
        if (i < _gridSize[0] - 1) {
          int neighSign = (_phi[ijkToid(i + 1, j, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i + 1, j, k)] - _phi[ijkToid(i, j, k)];
            double x = domainMin[0] - _h[0] * _phi[ijkToid(i, j, k)] / (theta) +
                       (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _ufaceLocation[faceijkToid(0, i + 1, j, k)] =
                Eigen::Array3d(x, y, z);
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
            _vfaceLocation[faceijkToid(1, i, j, k)] = Eigen::Array3d(x, y, z);
          }
        }
        if (j < _gridSize[1] - 1) {
          int neighSign = (_phi[ijkToid(i, j + 1, k)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i, j + 1, k)] - _phi[ijkToid(i, j, k)];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] - _h[1] * _phi[ijkToid(i, j, k)] / (theta) +
                       (j + 0.5) * _h[1];
            double z = domainMin[2] + (k + 0.5) * _h[2];
            _vfaceLocation[faceijkToid(1, i, j + 1, k)] =
                Eigen::Array3d(x, y, z);
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
            _wfaceLocation[faceijkToid(2, i, j, k)] = Eigen::Array3d(x, y, z);
          }
        }
        if (k < _gridSize[2] - 1) {
          int neighSign = (_phi[ijkToid(i, j, k + 1)] < 0) ? -1 : 1;
          if (centerSign != neighSign) {
            // Compute intersection location
            double theta = _phi[ijkToid(i, j, k + 1)] - _phi[ijkToid(i, j, k)];
            double x = domainMin[0] + (i + 0.5) * _h[0];
            double y = domainMin[1] + (j + 0.5) * _h[1];
            double z = domainMin[2] +
                       -_h[2] * _phi[ijkToid(i, j, k)] / (theta) +
                       (k + 0.5) * _h[2];
            _wfaceLocation[faceijkToid(2, i, j, k + 1)] =
                Eigen::Array3d(x, y, z);
          }
        }
      }
}

Eigen::Array3d DualCubes3::computeNormal(int cellId) {
  Eigen::Vector3d normal;
  auto h = getH();
  auto &phi = getScalarData(_phiId);
  double cell = phi[cellId];
  auto ijk = idToijk(cellId);
  int i = ijk[0], j = ijk[1], k = ijk[2];

  double neigh;
  if (i < _gridSize[0]) {
    neigh = phi[ijkToid(i + 1, j, k)];
    normal[0] = (neigh - cell) / h[0];
  } else {
    neigh = phi[ijkToid(i - 1, j, k)];
    normal[0] = (cell - neigh) / h[0];
  }

  if (j < _gridSize[1]) {
    neigh = phi[ijkToid(i, j + 1, k)];
    normal[1] = (neigh - cell) / h[1];
  } else {
    neigh = phi[ijkToid(i, j - 1, k)];
    normal[1] = (cell - neigh) / h[1];
  }

  if (k < _gridSize[2]) {
    neigh = phi[ijkToid(i, j, k + 1)];
    normal[2] = (neigh - cell) / h[2];
  } else {
    neigh = phi[ijkToid(i, j, k - 1)];
    normal[2] = (cell - neigh) / h[2];
  }
  return normal.normalized().array();
}

void DualCubes3::defineVelocity() {
  auto &_u = getFaceScalarData(0, _velocityId);
  auto &_v = getFaceScalarData(1, _velocityId);
  auto &_w = getFaceScalarData(2, _velocityId);
  // u velocities
  for (int id = 0; id < faceCount(0); id++) {
    auto p = facePosition(0, id);
    _u[id] = p[1];
  }

  // v velocities
  for (int id = 0; id < faceCount(1); id++) {
    auto p = facePosition(1, id);
    _v[id] = -p[0];
  }

  // w velocities
  for (int id = 0; id < faceCount(2); id++) {
    auto p = facePosition(2, id);
    _w[id] = 0;
  }
}

void DualCubes3::extractSurface() {
  Ramuh::DualMarching3 surface(_gridSize);
  auto &_phi = getScalarData("phi");
  auto &cellGradient = getArrayData("cellGradient");
  auto &_ufaceLocation = getFaceArrayData(0, "facePosition");
  auto &_vfaceLocation = getFaceArrayData(1, "facePosition");
  auto &_wfaceLocation = getFaceArrayData(2, "facePosition");
  auto &_ufaceNormals = getFaceArrayData(0, "faceNormal");
  auto &_vfaceNormals = getFaceArrayData(1, "faceNormal");
  auto &_wfaceNormals = getFaceArrayData(2, "faceNormal");

  std::vector<Ramuh::DualMarching3> surfaces(omp_get_max_threads(),
                                             Ramuh::DualMarching3(_gridSize));

  // #pragma omp parallel for
  for (int k = 0; k < _gridSize[2] - 1; k++) {
    for (int j = 0; j < _gridSize[1] - 1; j++) {
      for (int i = 0; i < _gridSize[0] - 1; i++) {
        std::vector<Eigen::Array3d> normalPosition;
        std::vector<Eigen::Vector3d> normal;
        bool isSurface = false;
        int id = ijkToid(i, j, k);

        // yz-plane
        if (signChange(_phi[ijkToid(i, j, k)], _phi[ijkToid(i + 1, j, k)])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[faceijkToid(0, i + 1, j, k)]);
          normalPosition.emplace_back(
              _ufaceLocation[faceijkToid(0, i + 1, j, k)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j + 1, k)],
                       _phi[ijkToid(i, j + 1, k)])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[faceijkToid(0, i + 1, j + 1, k)]);
          normalPosition.emplace_back(
              _ufaceLocation[faceijkToid(0, i + 1, j + 1, k)]);
        }
        if (signChange(_phi[ijkToid(i, j, k + 1)],
                       _phi[ijkToid(i + 1, j, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_ufaceNormals[faceijkToid(0, i + 1, j, k + 1)]);
          normalPosition.emplace_back(
              _ufaceLocation[faceijkToid(0, i + 1, j, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j + 1, k + 1)],
                       _phi[ijkToid(i, j + 1, k + 1)])) {
          isSurface = true;
          normal.emplace_back(
              _ufaceNormals[faceijkToid(0, i + 1, j + 1, k + 1)]);
          normalPosition.emplace_back(
              _ufaceLocation[faceijkToid(0, i + 1, j + 1, k + 1)]);
        }

        // xz-plane
        if (signChange(_phi[ijkToid(i + 1, j, k)],
                       _phi[ijkToid(i + 1, j + 1, k)])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[faceijkToid(1, i + 1, j + 1, k)]);
          normalPosition.emplace_back(
              _vfaceLocation[faceijkToid(1, i + 1, j + 1, k)]);
        }
        if (signChange(_phi[ijkToid(i, j + 1, k)], _phi[ijkToid(i, j, k)])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[faceijkToid(1, i, j + 1, k)]);
          normalPosition.emplace_back(
              _vfaceLocation[faceijkToid(1, i, j + 1, k)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j, k + 1)],
                       _phi[ijkToid(i + 1, j + 1, k + 1)])) {
          isSurface = true;
          normal.emplace_back(
              _vfaceNormals[faceijkToid(1, i + 1, j + 1, k + 1)]);
          normalPosition.emplace_back(
              _vfaceLocation[faceijkToid(1, i + 1, j + 1, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i, j + 1, k + 1)],
                       _phi[ijkToid(i, j, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_vfaceNormals[faceijkToid(1, i, j + 1, k + 1)]);
          normalPosition.emplace_back(
              _vfaceLocation[faceijkToid(1, i, j + 1, k + 1)]);
        }

        // xy-plane
        if (signChange(_phi[ijkToid(i, j, k)], _phi[ijkToid(i, j, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[faceijkToid(2, i, j, k + 1)]);
          normalPosition.emplace_back(
              _wfaceLocation[faceijkToid(2, i, j, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j, k)],
                       _phi[ijkToid(i + 1, j, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[faceijkToid(2, i + 1, j, k + 1)]);
          normalPosition.emplace_back(
              _wfaceLocation[faceijkToid(2, i + 1, j, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i + 1, j + 1, k)],
                       _phi[ijkToid(i + 1, j + 1, k + 1)])) {
          isSurface = true;
          normal.emplace_back(
              _wfaceNormals[faceijkToid(2, i + 1, j + 1, k + 1)]);
          normalPosition.emplace_back(
              _wfaceLocation[faceijkToid(2, i + 1, j + 1, k + 1)]);
        }
        if (signChange(_phi[ijkToid(i, j + 1, k)],
                       _phi[ijkToid(i, j + 1, k + 1)])) {
          isSurface = true;
          normal.emplace_back(_wfaceNormals[faceijkToid(2, i, j + 1, k + 1)]);
          normalPosition.emplace_back(
              _wfaceLocation[faceijkToid(2, i, j + 1, k + 1)]);
        }
        // Solve quadratic function
        if (isSurface) {
          Eigen::Array3d cubeMin =
              _domain.min() +
              _h.cwiseProduct(Eigen::Array3d(i + 0.5, j + 0.5, k + 0.5));
          Eigen::Array3d cubeMax =
              _domain.min() +
              _h.cwiseProduct(Eigen::Array3d(i + 1.5, j + 1.5, k + 1.5));
          int thread = omp_get_thread_num();
          auto x = surfaces[thread].evaluateCube(
              std::make_tuple(i, j, k), normalPosition, normal,
              Ramuh::BoundingBox3(cubeMin, cubeMax));
          // std::cout << x[0] << " " << x[1] << " " << x[2] <<
          // std::endl;
        }
      }
    }
  }
  for (size_t i = 0; i < omp_get_max_threads(); i++) {
    surface.merge(surfaces[i]);
  }

  surface.reconstruct();
}

bool DualCubes3::signChange(double valueA, double valueB) {
  int signA = (valueA < 0) ? -1 : 1;
  int signB = (valueB < 0) ? -1 : 1;
  return signA != signB;
}
} // namespace Carbuncle