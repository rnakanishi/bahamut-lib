#include <structures/grid3.h>
#include <geometry/matrix.h>
#include <utils/macros.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <queue>
#include <set>
#include <omp.h>

namespace Ramuh {

RegularGrid3::RegularGrid3() {
  setResolution(Vector3i(32, 32, 32));
  _h = Vector3d(1. / 32);
  _domainSize = Vector3d(1.0);
  _dt = 0.01;
}

Vector3d RegularGrid3::domainSize() { return _domainSize; }

Vector3i RegularGrid3::resolution() { return _resolution; }

Vector3d RegularGrid3::h() { return _h; }

Vector3d RegularGrid3::gridSize() { return _domainSize; }

void RegularGrid3::setSize(Vector3d newSize) {
  _domainSize = newSize;
  setH(_domainSize / _resolution);
}

void RegularGrid3::setResolution(Vector3i newResolution) {
  _resolution = newResolution;
  setH(_domainSize / _resolution);

  _u.resize(_resolution.x() + 1);
  for (auto &row : _u) {
    row.resize(_resolution.y());
    for (auto &depth : row)
      depth.resize(_resolution.z());
  }

  _v.resize(_resolution.x());
  for (auto &row : _v) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z());
  }

  _w.resize(_resolution.x());
  for (auto &row : _w) {
    row.resize(_resolution.y());
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
  }

  _material.resize(_resolution.x());
  for (auto &row : _material) {
    row.resize(_resolution.y());
    for (auto &depth : row)
      depth.resize(_resolution.z(), Material::FluidMaterial::AIR);
  }
}

void RegularGrid3::setH(Vector3d newH) { _h = newH; }

void RegularGrid3::macComarckVelocityAdvection() {
  Matrix3<Vector3d> utemp, vtemp, wtemp;
  utemp.changeSize(
      Vector3i(_resolution[0] + 1, _resolution[1], _resolution[2]));
  vtemp.changeSize(
      Vector3i(_resolution[0], _resolution[1] + 1, _resolution[2]));
  wtemp.changeSize(
      Vector3i(_resolution[0], _resolution[1], _resolution[2] + 1));

// Convective term for velocity U
#pragma omp parellel for
  for (int i = 1; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
        Eigen::Vector3d velocity;
        Eigen::Array3i index;
        double newU = _u[i][j][k][0];

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(0, h[1] / 2, h[2] / 2);
        velocity[0] = _u[i][j][k][0];
        velocity[1] = _interpolateVelocityV(position);
        velocity[2] = _interpolateVelocityW(position);
        position -= velocity.array() * _dt;
        index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();

        if ((velocity.norm() > 1e-8) && index[0] >= 0 && index[1] >= 0 &&
            index[2] >= 0 && index[0] < _resolution[0] &&
            index[1] < _resolution[1] && index[2] < _resolution[2]) {
          double u_n;
          u_n = _interpolateVelocityU(position);
          velocity[0] = _interpolateVelocityU(position);
          velocity[1] = _interpolateVelocityV(position);
          velocity[2] = _interpolateVelocityW(position);
          position += velocity.array() * _dt;
          double u_n1_hat = _interpolateVelocityU(position);
          double error = 0.5 * (_u[i][j][k][0] - u_n1_hat);
          newU = u_n + error;
        } else {
          newU = _u[i][j][k][0];
        }
        utemp[i][j][k] = newU;
      }

// Convective term for velocity V
#pragma omp parellel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 1; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
        Eigen::Vector3d velocity;
        Eigen::Array3i index;
        double newV = _v[i][j][k][1];

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(h[0] / 2, 0, h[2] / 2);
        velocity[0] = _interpolateVelocityU(position);
        velocity[1] = _v[i][j][k][1];
        velocity[2] = _interpolateVelocityW(position);
        position -= velocity.array() * _dt;
        index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();

        if ((velocity.norm() > 1e-8) && index[0] >= 0 && index[1] >= 0 &&
            index[2] >= 0 && index[0] < _resolution[0] &&
            index[1] < _resolution[1] && index[2] < _resolution[2]) {
          double v_n;
          v_n = _interpolateVelocityV(position);
          velocity[0] = _interpolateVelocityU(position);
          velocity[1] = _interpolateVelocityV(position);
          velocity[2] = _interpolateVelocityW(position);
          position += velocity.array() * _dt;
          double v_n1_hat = _interpolateVelocityV(position);
          double error = 0.5 * (_v[i][j][k][1] - v_n1_hat);
          newV = v_n + error;
        } else {
          newV = _v[i][j][k][1];
        }
        vtemp[i][j][k] = newV;
      }

// Convective term for velocity W
#pragma omp parellel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 1; k < _resolution.z(); k++) {
        Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
        Eigen::Vector3d velocity;
        Eigen::Array3i index;
        double newW = _w[i][j][k][2];

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(h[0] / 2, h[1] / 2, 0);
        velocity[0] = _interpolateVelocityU(position);
        velocity[1] = _interpolateVelocityV(position);
        velocity[2] = _w[i][j][k][2];
        position -= velocity.array() * _dt;
        index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();

        if ((velocity.norm() > 1e-8) && index[0] >= 0 && index[1] >= 0 &&
            index[2] >= 0 && index[0] < _resolution[0] &&
            index[1] < _resolution[1] && index[2] < _resolution[2]) {
          double w_n;
          w_n = _interpolateVelocityW(position);
          velocity[0] = _interpolateVelocityU(position);
          velocity[1] = _interpolateVelocityV(position);
          velocity[2] = _interpolateVelocityW(position);
          position += velocity.array() * _dt;
          double w_n1_hat = _interpolateVelocityW(position);
          double error = 0.5 * (_w[i][j][k][2] - w_n1_hat);
          newW = w_n + error;
        } else {
          newW = _w[i][j][k][2];
        }
        wtemp[i][j][k] = newW;
      }

#pragma omp parallel for
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _u[i][j][k][0] = utemp[i][j][k][0];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _v[i][j][k][1] = vtemp[i][j][k][1];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        _w[i][j][k][2] = wtemp[i][j][k][2];
      }
}

void RegularGrid3::advectGridVelocity() {
  // TODO: change to dynamic allocation
  Matrix3<Vector3d> utemp, vtemp, wtemp;

  utemp.changeSize(
      Vector3i(_resolution.x() + 1, _resolution.y(), _resolution.z()));
  vtemp.changeSize(
      Vector3i(_resolution.x(), _resolution.y() + 1, _resolution.z()));
  wtemp.changeSize(
      Vector3i(_resolution.x(), _resolution.y(), _resolution.z() + 1));

// // Interpolate other components
#pragma omp parallel for
  // x faces - all velocities components
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        double vVel = 0.0, wVel = 0.0;
        int vCount = 0, wCount = 0;

        if (i > 0) {
          vCount += 2;
          wCount += 2;
          vVel += _v[i - 1][j][k].y() + _v[i - 1][j + 1][k].y();
          wVel += _w[i - 1][j][k].z() + _w[i - 1][j][k + 1].z();
        }
        if (i < _resolution.x()) {
          vCount += 2;
          wCount += 2;
          vVel += _v[i][j][k].y() + _v[i][j + 1][k].y();
          wVel += _w[i][j][k].z() + _w[i][j][k + 1].z();
        }
        utemp[i][j][k].x(_u[i][j][k].x());
        utemp[i][j][k].y(vVel / vCount);
        utemp[i][j][k].z(wVel / wCount);
      }

#pragma omp parallel for
  // y faces - all velocities components
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        double uVel = 0.0, wVel = 0.0;
        int uCount = 0, wCount = 0;

        if (j > 0) {
          uCount += 2;
          wCount += 2;
          uVel += _u[i][j - 1][k].x() + _u[i + 1][j - 1][k].x();
          wVel += _w[i][j - 1][k].z() + _w[i][j - 1][k + 1].z();
        }
        if (j < _resolution.y()) {
          uCount += 2;
          wCount += 2;
          uVel += _u[i][j][k].x() + _u[i + 1][j][k].x();
          wVel += _w[i][j][k].z() + _w[i][j][k + 1].z();
        }
        vtemp[i][j][k].x(uVel / uCount);
        vtemp[i][j][k].y(_v[i][j][k].y());
        vtemp[i][j][k].z(wVel / wCount);
      }

#pragma omp parallel for
  // z faces - all velocities components
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        double uVel = 0.0, vVel = 0.0;
        int uCount = 0, vCount = 0;

        if (j > 0) {
          uCount += 2;
          vCount += 2;
          uVel += _u[i][j][k - 1].x() + _u[i + 1][j][k - 1].x();
          vVel += _v[i][j][k - 1].z() + _v[i][j + 1][k - 1].z();
        }
        if (j < _resolution.y()) {
          uCount += 2;
          vCount += 2;
          uVel += _u[i][j][k].x() + _u[i + 1][j][k].x();
          vVel += _v[i][j][k].z() + _v[i][j + 1][k].z();
        }
        wtemp[i][j][k].x(uVel / uCount);
        wtemp[i][j][k].y(vVel / vCount);
        wtemp[i][j][k].z(_w[i][j][k].z());
      }

// Solve convection term (material derivative) - u velocities
#pragma omp parallel for
  for (int i = 1; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Vector3d position, backPosition, velocity, cellCenter;
        Vector3d h(_h.x(), _h.y(), _h.z());
        Vector3i index;
        double newVelocity = _u[i][j][k].x();

        velocity = _u[i][j][k];
        position = Vector3d(i, j, k) * h + Vector3d(0, h.y() / 2, h.z() / 2);
        backPosition = position - velocity * _dt;
        index.x(std::floor(backPosition.x() / h.x()));
        index.y(std::floor(backPosition.y() / h.y()));
        index.z(std::floor(backPosition.z() / h.z()));

        if (velocity.length() > 1e-8 && index > Vector3i(0) &&
            index < _resolution - Vector3i(1)) {
          // back position fall inside domain
          cellCenter = Vector3d(index.x(), index.y(), index.z()) * h + h / 2.0;
          std::vector<int> iCand, jCand, kCand;
          iCand.push_back(index.x());
          iCand.push_back(index.x() + 1);
          jCand.push_back(index.y());
          kCand.push_back(index.z());
          if (backPosition.y() > cellCenter.y() &&
              index.y() < _resolution.y() - 1)
            jCand.push_back(index.y() + 1);
          else if (backPosition.y() < cellCenter.y() && index.y() > 0)
            jCand.push_back(index.y() - 1);
          if (backPosition.z() > cellCenter.z() &&
              index.z() < _resolution.z() - 1)
            kCand.push_back(index.z() + 1);
          else if (backPosition.z() < cellCenter.z() && index.z() > 0)
            kCand.push_back(index.z() - 1);

          newVelocity = 0.0;
          double distanceCount = 0.0, distance = 0.0;
          for (auto u : iCand)
            for (auto v : jCand)
              for (auto w : kCand) {
                position =
                    Vector3d(u, v, w) * h + Vector3d(0, h.y() / 2, h.z() / 2);
                distance = (backPosition - position).length();
                distanceCount += 1. / distance;
                newVelocity += (utemp[u][v][w].x() / distance);
              }
          newVelocity /= distanceCount;
        }
        _u[i][j][k].x(newVelocity);
      }

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 1; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Vector3d position, backPosition, velocity, cellCenter;
        Vector3d h(_h.x(), _h.y(), _h.z());
        Vector3i index;
        double newVelocity = _v[i][j][k].y();

        velocity = _v[i][j][k];
        position = Vector3d(i, j, k) * h + Vector3d(h.x() / 2, 0, h.z() / 2);
        backPosition = position - velocity * _dt;
        index.x(std::floor(backPosition.x() / h.x()));
        index.y(std::floor(backPosition.y() / h.y()));
        index.z(std::floor(backPosition.z() / h.z()));

        if (velocity.length() > 1e-8 && index > Vector3i(0) &&
            index < _resolution - Vector3i(1)) {
          // back position fall inside domain
          cellCenter = Vector3d(index.x(), index.y(), index.z()) * h + h / 2.0;
          std::vector<int> iCand, jCand, kCand;
          iCand.push_back(index.x());
          jCand.push_back(index.y());
          jCand.push_back(index.y() + 1);
          kCand.push_back(index.z());
          if (backPosition.x() > cellCenter.x() &&
              index.x() < _resolution.x() - 1)
            iCand.push_back(index.x() + 1);
          else if (backPosition.x() < cellCenter.x() && index.x() > 0)
            iCand.push_back(index.x() - 1);
          if (backPosition.z() > cellCenter.z() &&
              index.z() < _resolution.z() - 1)
            kCand.push_back(index.z() + 1);
          else if (backPosition.z() < cellCenter.z() && index.z() > 0)
            kCand.push_back(index.z() - 1);

          newVelocity = 0.0;
          double distanceCount = 0.0, distance = 0.0;
          for (auto u : iCand)
            for (auto v : jCand)
              for (auto w : kCand) {
                position =
                    Vector3i(u, v, w) * h + Vector3d(h.x() / 2, 0, h.z() / 2);
                distance = (backPosition - position).length();
                distanceCount += 1. / distance;
                newVelocity += (vtemp[u][v][w].y() / distance);
              }
          newVelocity /= distanceCount;
        }
        _v[i][j][k].y(newVelocity);
      }

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 1; k < _resolution.z() + 1; k++) {
        Vector3d position, backPosition, velocity, cellCenter;
        Vector3d h(_h.x(), _h.y(), _h.z());
        Vector3i index;
        double newVelocity = _w[i][j][k].z();

        velocity = _w[i][j][k];
        position = Vector3d(i, j, k) * h + Vector3d(h.x(), h.y(), 0) / 2;
        backPosition = position - velocity * _dt;
        index.x(std::floor(backPosition.x() / h.x()));
        index.y(std::floor(backPosition.y() / h.y()));
        index.z(std::floor(backPosition.z() / h.z()));

        if (velocity.length() > 1e-8 && index > Vector3i(0) &&
            index < _resolution - Vector3i(1)) {
          cellCenter = Vector3d(index.x(), index.y(), index.z()) * h + h / 2.0;
          std::vector<int> iCand, jCand, kCand;
          iCand.push_back(index.x());
          jCand.push_back(index.y());
          kCand.push_back(index.z());
          kCand.push_back(index.z() + 1);
          if (backPosition.x() > cellCenter.x() &&
              index.x() < _resolution.x() - 1)
            iCand.push_back(index.x() + 1);
          else if (backPosition.x() < cellCenter.x() && index.x() > 0)
            iCand.push_back(index.x() - 1);
          if (backPosition.y() > cellCenter.y() &&
              index.y() < _resolution.y() - 1)
            jCand.push_back(index.y() + 1);
          else if (backPosition.y() < cellCenter.y() && index.y() > 0)
            jCand.push_back(index.y() - 1);

          newVelocity = 0.0;
          double distanceCount = 0.0, distance = 0.0;
          for (auto u : iCand)
            for (auto v : jCand)
              for (auto w : kCand) {
                position =
                    Vector3i(u, v, w) * h + Vector3d(h.x(), h.y(), 0) / 2;
                distance = (backPosition - position).length();
                distanceCount += 1. / distance;
                newVelocity += (wtemp[u][v][w].z() / distance);
              }
          newVelocity /= distanceCount;
        }
        _w[i][j][k].z(newVelocity);
      }
}

void RegularGrid3::setVelocity() {

  for (int i = 0; i < _resolution.x() + 1; i++) {
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _u[i][j][k] = Vector3d(0, 0, 0);
      }
  }

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _v[i][j][k] = Vector3d(0, 0, 0);
      }
  }

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        _w[i][j][k] = Vector3d(0, 0, 0);
      }
  }
}

int RegularGrid3::cellCount() {
  return _resolution.x() * _resolution.y() * _resolution.z();
}

int RegularGrid3::ijkToId(int i, int j, int k) {
  return k * _resolution.x() * _resolution.y() + j * _resolution.x() + i;
}

Vector3i RegularGrid3::idToijk(int id) { NOT_IMPLEMENTED(); }

void RegularGrid3::printFaceVelocity() {

  // std::cerr << "==== u: \n";
  // for (int k = 0; k < _resolution.z(); k++) {
  //   for (int j = 0; j < _resolution.y(); j++) {
  //     for (int i = 0; i < _resolution.x() + 1; i++) {
  //       std::cerr << _u[i][j][k].x() << " ";
  //     }
  //     std::cerr << std::endl;
  //   }
  //   std::cerr << std::endl;
  // }
  std::cerr << "==== v: \n";
  int k = _resolution.y() / 2;
  // for (int k = 0; k < _resolution.z(); k++) {
  for (int j = 0; j < _resolution.y() + 1; j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      std::cerr << _v[i][j][k].y() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  // }
  // std::cerr << "==== w: \n";
  // for (int k = 0; k < _resolution.z(); k++) {
  //   for (int j = 0; j < _resolution.y() + 1; j++) {
  //     for (int i = 0; i < _resolution.x(); i++) {
  //       std::cerr << _w[i][j][k].z() << " ";
  //     }
  //     std::cerr << std::endl;
  //   }
  //   std::cerr << std::endl;
  // }
}

void RegularGrid3::boundaryVelocities() {
  // Z velocities
  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      _w[i][j][0].z(0.0);
      _w[i][j][_resolution.z()].z(0.0);
    }
  }
  for (int k = 0; k < _resolution.z(); k++) {
    for (int i = 0; i < _resolution.x(); i++) {
      _v[i][0][k].y(0.0);
      _v[i][_resolution.y()][k].y(0.0);
    }
  }
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      _u[0][j][k].x(0.0);
      _u[_resolution.x()][j][k].x(0.0);
    }
  }
}

void RegularGrid3::addGravity() {
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 1; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        // if (_material[i][j - 1][k] == Material::FluidMaterial::FLUID) {
        double vel = _v[i][j][k].y() - 9.81 * _dt;
        _v[i][j][k].y(vel);
      }
    }
  }
}

void RegularGrid3::extrapolateVelocity() {
  std::queue<Vector3i> processingCells;
  // TODO: Change processed cells to check phi instead of construct distance vec
  Matrix3<int> processedCells;
  processedCells.changeSize(_resolution);

  // Extrapolate velocity from fluid to air cells

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int k = 0; k < _resolution.z(); k++) {
        // Check if is a air surface cell
        if (_material[i][j][k] == Material::FluidMaterial::FLUID) {
          processedCells[i][j][k] = 0;
          // Look for air neighbors
          // Add air cells to processing queue
          if (i > 0 && _material[i - 1][j][k] == Material::FluidMaterial::AIR) {
            processedCells[i - 1][j][k] = 1;
            processingCells.push(Vector3i(i - 1, j, k));
            _u[i - 1][j][k].x(_u[i][j][k].x());
            _v[i - 1][j][k].y(_v[i][j][k].y());
            _v[i - 1][j + 1][k].y(_v[i][j + 1][k].y());
            _w[i - 1][j][k].z(_w[i][j][k].z());
            _w[i - 1][j][k + 1].z(_w[i][j][k + 1].z());
          }
          if (i < _resolution.x() - 1 &&
              _material[i + 1][j][k] == Material::FluidMaterial::AIR) {
            processedCells[i + 1][j][k] = 1;
            processingCells.push(Vector3i(i + 1, j, k));
            _u[i + 2][j][k].x(_u[i + 1][j][k].x());
            _v[i + 1][j][k].y(_v[i][j][k].y());
            _v[i + 1][j + 1][k].y(_v[i][j + 1][k].y());
            _w[i + 1][j][k].z(_w[i][j][k].z());
            _w[i + 1][j][k + 1].z(_w[i][j][k + 1].z());
          }
          if (j > 0 && _material[i][j - 1][k] == Material::FluidMaterial::AIR) {
            processedCells[i][j - 1][k] = 1;
            processingCells.push(Vector3i(i, j - 1, k));
            _v[i][j - 1][k].y(_v[i][j][k].y());
            _u[i][j - 1][k].x(_u[i][j][k].x());
            _u[i + 1][j - 1][k].x(_u[i + 1][j][k].x());
            _w[i][j - 1][k].z(_w[i][j][k].z());
            _w[i][j - 1][k + 1].z(_w[i][j][k + 1].z());
          }
          if (j < _resolution.y() - 1 &&
              _material[i][j + 1][k] == Material::FluidMaterial::AIR) {
            processedCells[i][j + 1][k] = 1;
            processingCells.push(Vector3i(i, j + 1, k));
            _v[i][j + 2][k].y(_v[i][j + 1][k].y());
            _u[i][j + 1][k].x(_u[i][j][k].x());
            _u[i + 1][j + 1][k].x(_u[i + 1][j][k].x());
            _w[i][j + 1][k].z(_w[i][j][k].z());
            _w[i][j + 1][k + 1].z(_w[i][j][k + 1].z());
          }
          if (k > 0 && _material[i][j][k - 1] == Material::FluidMaterial::AIR) {
            processedCells[i][j][k - 1] = 1;
            processingCells.push(Vector3i(i, j, k - 1));
            _w[i][j][k - 1].z(_w[i][j][k].z());
            _u[i][j][k - 1].x(_u[i][j][k].x());
            _u[i + 1][j][k - 1].x(_u[i + 1][j][k].x());
            _v[i][j][k - 1].y(_v[i][j][k].y());
            _v[i][j + 1][k - 1].y(_v[i][j + 1][k].y());
          }
          if (k < _resolution.y() - 1 &&
              _material[i][j][k + 1] == Material::FluidMaterial::AIR) {
            processedCells[i][j][k + 1] = 1;
            processingCells.push(Vector3i(i, j, k + 1));
            _w[i][j][k + 2].z(_w[i][j][k + 1].z());
            _u[i][j][k + 1].x(_u[i][j][k].x());
            _u[i + 1][j][k + 1].x(_u[i + 1][j][k].x());
            _v[i][j][k + 1].y(_v[i][j][k].y());
            _v[i][j + 1][k + 1].y(_v[i][j + 1][k].y());
          }
        }
      }
    }
  }
  while (!processingCells.empty()) {
    Vector3i currCell = processingCells.front();
    processingCells.pop();
    int i = currCell.x(), j = currCell.y(), k = currCell.z();

    // Check for neighbor airs
    if (i > 0 && _material[i - 1][j][k] == Material::FluidMaterial::AIR) {
      if (processedCells[i - 1][j][k] > processedCells[i][j][k] + 1) {
        processedCells[i - 1][j][k] = processedCells[i][j][k] + 1;
        processingCells.push(Vector3i(i - 1, j, 0));
        _u[i - 1][j][k].x(_u[i][j][k].x());
        _v[i - 1][j][k].y(_v[i][j][k].y());
        _v[i - 1][j + 1][k].y(_v[i][j + 1][k].y());
        _w[i - 1][j][k].z(_w[i][j][k].z());
        _w[i - 1][j][k + 1].z(_w[i][j][k + 1].z());
      }
    }
    if (i < _resolution.x() - 1 &&
        _material[i + 1][j][k] == Material::FluidMaterial::AIR) {
      if (processedCells[i + 1][j][k] > processedCells[i][j][k] + 1) {
        processedCells[i + 1][j][k] = processedCells[i][j][k] + 1;
        processingCells.push(Vector3i(i + 1, j, 0));
        _u[i + 2][j][k].x(_u[i + 1][j][k].x());
        _v[i + 1][j][k].y(_v[i][j][k].y());
        _v[i + 1][j + 1][k].y(_v[i][j + 1][k].y());
        _w[i + 1][j][k].z(_w[i][j][k].z());
        _w[i + 1][j][k + 1].z(_w[i][j][k + 1].z());
      }
    }
    if (j > 0 && _material[i][j - 1][k] == Material::FluidMaterial::AIR) {
      if (processedCells[i][j - 1][k] > processedCells[i][j][k] + 1) {
        processedCells[i][j - 1][k] = processedCells[i][j][k] + 1;
        processingCells.push(Vector3i(i, j - 1, 0));
        _v[i][j - 1][k].y(_v[i][j][k].y());
        _u[i][j - 1][k].x(_u[i][j][k].x());
        _u[i + 1][j - 1][k].x(_u[i + 1][j][k].x());
        _w[i][j - 1][k].z(_w[i][j][k].z());
        _w[i][j - 1][k + 1].z(_w[i][j][k + 1].z());
      }
    }
    if (j < _resolution.y() - 1 &&
        _material[i][j + 1][k] == Material::FluidMaterial::AIR) {
      if (processedCells[i][j + 1][k] > processedCells[i][j][k] + 1) {
        processedCells[i][j + 1][k] = processedCells[i][j][k] + 1;
        processingCells.push(Vector3i(i, j + 1, 0));
        _v[i][j + 2][k].y(_v[i][j + 1][k].y());
        _u[i][j + 1][k].x(_u[i][j][k].x());
        _u[i + 1][j + 1][k].x(_u[i + 1][j][k].x());
        _w[i][j + 1][k].z(_w[i][j][k].z());
        _w[i][j + 1][k + 1].z(_w[i][j][k + 1].z());
      }
      if (k > 0 && _material[i][j][k - 1] == Material::FluidMaterial::AIR) {
        if (processedCells[i][j][k + 1] > processedCells[i][j][k] + 1) {
          processedCells[i][j][k - 1] = 1;
          processingCells.push(Vector3i(i, j, k - 1));
          _w[i][j][k - 1].z(_w[i][j][k].z());
          _u[i][j][k - 1].x(_u[i][j][k].x());
          _u[i + 1][j][k - 1].x(_u[i + 1][j][k].x());
          _v[i][j][k - 1].y(_v[i][j][k].y());
          _v[i][j + 1][k - 1].y(_v[i][j + 1][k].y());
        }
      }
      if (k < _resolution.y() - 1 &&
          _material[i][j][k + 1] == Material::FluidMaterial::AIR) {
        if (processedCells[i][j][k + 1] > processedCells[i][j][k] + 1) {
          processedCells[i][j][k + 1] = 1;
          processingCells.push(Vector3i(i, j, k + 1));
          _w[i][j][k + 2].z(_w[i][j][k + 1].z());
          _u[i][j][k + 1].x(_u[i][j][k].x());
          _u[i + 1][j][k + 1].x(_u[i + 1][j][k].x());
          _v[i][j][k + 1].y(_v[i][j][k].y());
          _v[i][j + 1][k + 1].y(_v[i][j + 1][k].y());
        }
      }
    }
  }
}

void RegularGrid3::solvePressure() {

  // Compute velocity divergent over cell center
  // std::vector<double> divergent(cellCount(), 0.0);
  int nCells = cellCount();
  Eigen::VectorXd divergent, pressure;
  Eigen::SparseMatrix<double> pressureMatrix(nCells + 1, nCells + 1);
  std::vector<Eigen::Triplet<double>> triplets;

  divergent = Eigen::VectorXd::Zero(nCells + 1);
  pressure = Eigen::VectorXd::Zero(nCells + 1);

// Solve pressure Poisson equation

#pragma omp parallel for
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        if (_material[i][j][k] != Material::FluidMaterial::FLUID)
          continue;

        int validCells = 0;
        int id = ijkToId(i, j, k);
        std::vector<Eigen::Triplet<double>> threadTriplet;

        if (i > 0) {
          validCells++;
          if (_material[i - 1][j][k] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijkToId(i - 1, j, k),
                                       -1 / (_h.x() * _h.x()));
        }
        if (i < _resolution.x() - 1) {
          validCells++;
          if (_material[i + 1][j][k] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijkToId(i + 1, j, k),
                                       -1 / (_h.x() * _h.x()));
        }
        if (j > 0) {
          validCells++;
          if (_material[i][j - 1][k] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijkToId(i, j - 1, k),
                                       -1 / (_h.x() * _h.x()));
        }
        if (j < _resolution.y() - 1) {
          validCells++;
          if (_material[i][j + 1][k] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijkToId(i, j + 1, k),
                                       -1 / (_h.x() * _h.x()));
        }
        if (k > 0) {
          validCells++;
          if (_material[i][j][k - 1] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijkToId(i, j, k - 1),
                                       -1 / (_h.x() * _h.x()));
        }
        if (k < _resolution.z() - 1) {
          validCells++;
          if (_material[i][j][k + 1] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijkToId(i, j, k + 1),
                                       -1 / (_h.x() * _h.x()));
        }

        threadTriplet.emplace_back(id, id, validCells / (_h.x() * _h.x()));
#pragma omp critical
        {
          triplets.insert(triplets.end(), threadTriplet.begin(),
                          threadTriplet.end());
        }

        divergent[id] = 0;
        divergent[id] -= (_u[i + 1][j][k].x() - _u[i][j][k].x()) / _h.x();
        divergent[id] -= (_v[i][j + 1][k].y() - _v[i][j][k].y()) / _h.y();
        divergent[id] -= (_w[i][j][k + 1].z() - _w[i][j][k].z()) / _h.z();
        divergent[id] /= _dt;
      }
    }
  }
  triplets.emplace_back(nCells, nCells, 1);

  pressureMatrix.setFromTriplets(triplets.begin(), triplets.end());

  // SOlve pressure Poisson system
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;
  solver.compute(pressureMatrix);
  pressure = solver.solve(divergent);

// Correct velocity through pressure gradient
#pragma omp parallel for
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 1; i < _resolution.x(); i++) {
        _u[i][j][k].x(_u[i][j][k].x() -
                      _dt * (pressure[ijkToId(i, j, k)] -
                             pressure[ijkToId(i - 1, j, k)]) /
                          _h.x());
      }
    }
  }
#pragma omp parallel for
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 1; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        _v[i][j][k].y(_v[i][j][k].y() -
                      _dt * (pressure[ijkToId(i, j, k)] -
                             pressure[ijkToId(i, j - 1, k)]) /
                          _h.y());
      }
    }
  }
#pragma omp parallel for
  for (int k = 1; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        _w[i][j][k].z(_w[i][j][k].z() -
                      _dt * (pressure[ijkToId(i, j, k)] -
                             pressure[ijkToId(i, j, k - 1)]) /
                          _h.z());
      }
    }
  }
}

double RegularGrid3::_interpolateVelocityU(Eigen::Array3d position) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution.x(), _resolution.y(), _resolution.z());
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  std::vector<int> iCandidates, jCandidates, kCandidates;
  iCandidates.push_back(index[0]);
  iCandidates.push_back(index[0] + 1);
  if (index[1] < resolution[1])
    jCandidates.push_back(index[1]);
  if (index[2] < resolution[2])
    kCandidates.push_back(index[2]);
  if (position[1] > cellCenter[1] && index[1] < _resolution[1] - 1)
    jCandidates.push_back(index[1] + 1);
  else if (position[1] < cellCenter[1] && index[1] > 0)
    jCandidates.push_back(index[1] - 1);
  if (position[2] > cellCenter[2] && index[2] < _resolution[2] - 1)
    kCandidates.push_back(index[2] + 1);
  else if (position[2] < cellCenter[2] && index[2] > 0)
    kCandidates.push_back(index[2] - 1);

  double velocity(0);
  // Catmull-Rom like interpolation of u
  double distance = 0., distanceCount = 0., newVelocity = 1e8;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        Eigen::Array3d centerPosition =
            Eigen::Array3d(u, v, w) * h + Eigen::Array3d(0, h[1] / 2, h[2] / 2);
        if ((position - centerPosition).matrix().norm() < 1e-7)
          return _u[u][v][w].x();
        distance = (position - centerPosition).matrix().norm();
        distanceCount += 1. / distance;
        velocity += _u[u][v][w].x() / distance;
      }
  return velocity / distanceCount;
}

double RegularGrid3::_interpolateVelocityV(Eigen::Array3d position) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution.x(), _resolution.y(), _resolution.z());
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  std::vector<int> iCandidates, jCandidates, kCandidates;
  if (index[0] < resolution[0])
    iCandidates.push_back(index[0]);
  jCandidates.push_back(index[1]);
  jCandidates.push_back(index[1] + 1);
  if (index[2] < resolution[2])
    kCandidates.push_back(index[2]);
  if (position[0] > cellCenter[0] && index[0] < resolution[0] - 1)
    iCandidates.push_back(index[0] + 1);
  else if (position[0] < cellCenter[0] && index[0] > 0)
    iCandidates.push_back(index[0] - 1);
  if (position[2] > cellCenter[2] && index[2] < _resolution[2] - 1)
    kCandidates.push_back(index[2] + 1);
  else if (position[2] < cellCenter[2] && index[2] > 0)
    kCandidates.push_back(index[2] - 1);

  double velocity(0);
  // Catmull-Rom like interpolation of v
  double distance = 0., distanceCount = 0., newVelocity = 1e8;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        Eigen::Array3d centerPosition =
            Eigen::Array3d(u, v, w) * h + Eigen::Array3d(h[0] / 2, 0, h[2] / 2);
        if ((position - centerPosition).matrix().norm() < 1e-7)
          return _v[u][v][w].y();
        distance = (position - centerPosition).matrix().norm();
        distanceCount += 1. / distance;
        velocity += _v[u][v][w].y() / distance;
      }
  return velocity / distanceCount;
}

double RegularGrid3::_interpolateVelocityW(Eigen::Array3d position) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution.x(), _resolution.y(), _resolution.z());
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  std::vector<int> iCandidates, jCandidates, kCandidates;
  if (index[0] < resolution[0])
    iCandidates.push_back(index[0]);
  if (index[1] < resolution[1])
    jCandidates.push_back(index[1]);
  kCandidates.push_back(index[2]);
  kCandidates.push_back(index[2] + 1);
  if (position[0] > cellCenter[0] && index[0] < resolution[0] - 1)
    iCandidates.push_back(index[0] + 1);
  else if (position[0] < cellCenter[0] && index[0] > 0)
    iCandidates.push_back(index[0] - 1);
  if (position[1] > cellCenter[1] && index[1] < _resolution[1] - 1)
    jCandidates.push_back(index[1] + 1);
  else if (position[1] < cellCenter[1] && index[1] > 0)
    jCandidates.push_back(index[1] - 1);

  double velocity(0);
  // Catmull-Rom like interpolation of u
  double distance = 0., distanceCount = 0., newVelocity = 1e8;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        Eigen::Array3d centerPosition =
            Eigen::Array3d(u, v, w) * h + Eigen::Array3d(0, h[1] / 2, h[2] / 2);
        if ((position - centerPosition).matrix().norm() < 1e-7)
          return _w[u][v][w].z();
        distance = (position - centerPosition).matrix().norm();
        distanceCount += 1. / distance;
        velocity += _w[u][v][w].z() / distance;
      }
  return velocity / distanceCount;
}

} // namespace Ramuh
