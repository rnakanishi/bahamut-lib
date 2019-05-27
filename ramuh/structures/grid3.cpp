#include <structures/grid3.h>
#include <blas/interpolator.h>
#include <geometry/matrix.h>
#include <utils/macros.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <omp.h>
#include <utils/timer.hpp>

namespace Ramuh {

RegularGrid3::RegularGrid3() {
  setResolution(Vector3i(32, 32, 32));
  _h = Vector3d(1. / 32);
  _domainSize = Vector3d(1.0);
  _maxVelocity[0] = _maxVelocity[1] = _maxVelocity[2] = -1e8;

  _ellapsedDt = 0.0;
  _originalDt = _dt = 1.0 / 60;
  _currBuffer = 0;
  _tolerance = 1e-12;
}

double RegularGrid3::tolerance() { return _tolerance; }

void RegularGrid3::setTolerance(double tol) { _tolerance = tol; }

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
  for (int buffer = 0; buffer < 2; buffer++) {
    _u[buffer].resize(_resolution.x() + 1);
    _uFaceMaterial.resize(_resolution.x() + 1);
    for (int i = 0; i < _resolution.x() + 1; i++) {
      _u[buffer][i].resize(_resolution.y());
      _uFaceMaterial[i].resize(_resolution.y());
      for (int j = 0; j < _resolution.y(); j++) {
        _uFaceMaterial[i][j].resize(_resolution.z());
        _u[buffer][i][j].resize(_resolution.z());
      }
    }

    _v[buffer].resize(_resolution.x());
    _vFaceMaterial.resize(_resolution.x());
    for (int i = 0; i < _resolution.x(); i++) {
      _vFaceMaterial[i].resize(_resolution.y() + 1);
      _v[buffer][i].resize(_resolution.y() + 1);
      for (int j = 0; j < _resolution.y() + 1; j++) {
        _v[buffer][i][j].resize(_resolution.z());
        _vFaceMaterial[i][j].resize(_resolution.z());
      }
    }

    _w[buffer].resize(_resolution.x());
    _wFaceMaterial.resize(_resolution.x());
    for (int i = 0; i < _resolution.y(); i++) {
      _w[buffer][i].resize(_resolution.y());
      _wFaceMaterial[i].resize(_resolution.y());
      for (int j = 0; j < _resolution.y(); j++) {
        _w[buffer][i][j].resize(_resolution.z() + 1);
        _wFaceMaterial[i][j].resize(_resolution.z() + 1);
      }
    }

    _material.resize(_resolution.x());
    for (auto &row : _material) {
      row.resize(_resolution.y());
      for (auto &depth : row)
        depth.resize(_resolution.z(), Material::FluidMaterial::AIR);
    }
  }
}

void RegularGrid3::setH(Vector3d newH) { _h = newH; }

void RegularGrid3::swapBuffers() { _currBuffer = (_currBuffer + 1) % 2; }

bool RegularGrid3::advanceTime() {
  _ellapsedDt += _dt;
  if (_ellapsedDt < _originalDt)
    return false;
  _ellapsedDt = 0.0;
  _dt = _originalDt;
  return true;
}

void RegularGrid3::cfl() {
  // Find biggest velocity

  Eigen::Vector3d maxVel = Eigen::Vector3d(0.0, 0.0, 0.0);
  for (int id = 0; id < cellCount(); id++) {
    Eigen::Array3i ijk = idToijk(id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    Eigen::Vector3d vel;
    vel[0] = (_u[_currBuffer][i][j][k] + _u[_currBuffer][i + 1][j][k]) / 2.0;
    vel[1] = (_v[_currBuffer][i][j][k] + _v[_currBuffer][i][j + 1][k]) / 2.0;
    vel[2] = (_w[_currBuffer][i][j][k] + _w[_currBuffer][i][j][k + 1]) / 2.0;

    if (maxVel.norm() < vel.norm() * 1.05)
      maxVel = vel;
  }

  // Check if cfl condition applies
  // Half timestep if so
  if (maxVel.norm() * (_originalDt - _ellapsedDt) > 3 * _h[0]) {
    _dt = _dt / 2;
  }
}

void RegularGrid3::macComarckVelocityAdvection() {
  Matrix3<double> utemp, vtemp, wtemp;
  utemp.changeSize(
      Vector3i(_resolution[0] + 1, _resolution[1], _resolution[2]));
  vtemp.changeSize(
      Vector3i(_resolution[0], _resolution[1] + 1, _resolution[2]));
  wtemp.changeSize(
      Vector3i(_resolution[0], _resolution[1], _resolution[2] + 1));

  double clamp[2]; // 0: min, 1: max
// Convective term for velocity U
#pragma omp parellel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
        Eigen::Vector3d velocity;
        Eigen::Array3i index;

        if (_uFaceMaterial[i][j][k] != Material::FluidMaterial::FLUID) {
          utemp[i][j][k] = _u[_currBuffer][i][j][k];
        } else {
          double newU = _u[_currBuffer][i][j][k];

          position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                     Eigen::Array3d(0, h[1] / 2, h[2] / 2);
          velocity[0] = _u[_currBuffer][i][j][k];
          velocity[1] = _interpolateVelocityV(position);
          velocity[2] = _interpolateVelocityW(position);
          position -= velocity.array() * _dt;
          index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();

          if ((velocity.norm() > _tolerance) && index[0] >= 0 &&
              index[1] >= 0 && index[2] >= 0 && index[0] < _resolution[0] &&
              index[1] < _resolution[1] && index[2] < _resolution[2]) {
            double u_n;
            u_n = _interpolateVelocityU(position, clamp[0], clamp[1]);
            velocity[0] = _interpolateVelocityU(position);
            velocity[1] = _interpolateVelocityV(position);
            velocity[2] = _interpolateVelocityW(position);
            position += velocity.array() * _dt;
            double u_n1_hat = _interpolateVelocityU(position);
            double error = 0.5 * (_u[_currBuffer][i][j][k] - u_n1_hat);
            newU = u_n + error;
          } else {
            newU = _u[_currBuffer][i][j][k];
          }
          utemp[i][j][k] = std::max(clamp[0], std::min(clamp[1], newU));
        }
        // utemp[i][j][k] = newU;
        //       }

        // // Convective term for velocity V
        if (_vFaceMaterial[i][j][k] != Material::FluidMaterial::FLUID) {
          vtemp[i][j][k] = _v[_currBuffer][i][j][k];
        } else {
          double newV = _v[_currBuffer][i][j][k];

          position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                     Eigen::Array3d(h[0] / 2, 0, h[2] / 2);
          velocity[0] = _interpolateVelocityU(position);
          velocity[1] = _v[_currBuffer][i][j][k];
          velocity[2] = _interpolateVelocityW(position);
          position -= velocity.array() * _dt;
          index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();

          if ((velocity.norm() > _tolerance) && index[0] >= 0 &&
              index[1] >= 0 && index[2] >= 0 && index[0] < _resolution[0] &&
              index[1] < _resolution[1] && index[2] < _resolution[2]) {
            double v_n;
            v_n = _interpolateVelocityV(position, clamp[0], clamp[1]);
            velocity[0] = _interpolateVelocityU(position);
            velocity[1] = _interpolateVelocityV(position);
            velocity[2] = _interpolateVelocityW(position);
            position += velocity.array() * _dt;
            double v_n1_hat = _interpolateVelocityV(position);
            double error = 0.5 * (_v[_currBuffer][i][j][k] - v_n1_hat);
            newV = v_n + error;
          } else {
            newV = _v[_currBuffer][i][j][k];
          }
          vtemp[i][j][k] = std::max(clamp[0], std::min(clamp[1], newV));
        }
        // vtemp[i][j][k] = newV;
        //       }

        // // Convective term for velocity W
        if (_wFaceMaterial[i][j][k] != Material::FluidMaterial::FLUID) {
          wtemp[i][j][k] = _w[_currBuffer][i][j][k];
        } else {
          double newW = _w[_currBuffer][i][j][k];

          position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                     Eigen::Array3d(h[0] / 2, h[1] / 2, 0);
          velocity[0] = _interpolateVelocityU(position);
          velocity[1] = _interpolateVelocityV(position);
          velocity[2] = _w[_currBuffer][i][j][k];
          position -= velocity.array() * _dt;
          index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();

          if ((velocity.norm() > _tolerance) && index[0] >= 0 &&
              index[1] >= 0 && index[2] >= 0 && index[0] < _resolution[0] &&
              index[1] < _resolution[1] && index[2] < _resolution[2]) {
            double w_n;
            w_n = _interpolateVelocityW(position, clamp[0], clamp[1]);
            velocity[0] = _interpolateVelocityU(position);
            velocity[1] = _interpolateVelocityV(position);
            velocity[2] = _interpolateVelocityW(position);
            position += velocity.array() * _dt;
            double w_n1_hat = _interpolateVelocityW(position);
            double error = 0.5 * (_w[_currBuffer][i][j][k] - w_n1_hat);
            newW = w_n + error;
          } else {
            newW = _w[_currBuffer][i][j][k];
          }
          wtemp[i][j][k] = std::max(clamp[0], std::min(clamp[1], newW));
          // wtemp[i][j][k] = newW;
        }
      }

#pragma omp barrier
#pragma omp parallel for
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _u[_currBuffer][i][j][k] = utemp[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _v[_currBuffer][i][j][k] = vtemp[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        _w[_currBuffer][i][j][k] = wtemp[i][j][k];
      }
}

void RegularGrid3::advectGridVelocity() {
  // TODO: change to dynamic allocation
  Matrix3<double> utemp, vtemp, wtemp;
  Eigen::Array3d h(_h[0], _h[1], _h[2]);

  utemp.changeSize(
      Vector3i(_resolution.x() + 1, _resolution.y(), _resolution.z()), -1e8);
  vtemp.changeSize(
      Vector3i(_resolution.x(), _resolution.y() + 1, _resolution.z()), -1e8);
  wtemp.changeSize(
      Vector3i(_resolution.x(), _resolution.y(), _resolution.z() + 1), -1e8);

  std::vector<std::tuple<int, int, int>> fluidFaces;
  for (auto id : _fluidCells) {
    Eigen::Array3i ijk = idToijk(id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];
    fluidFaces.emplace_back(std::make_tuple(i, j, k));
    if (i < _resolution[0] - 1 &&
        _material[i + 1][j][k] == Material::FluidMaterial::AIR) {
      fluidFaces.emplace_back(std::make_tuple(i + 1, j, k));
    }
  }
  // -------------------------------------------
  // Solve convection term (material derivative) - u velocities
  Eigen::Array3i ijk;

#pragma omp parallel for
  for (int it = 0; it < fluidFaces.size(); it++) {
    int i, j, k;
    Eigen::Array3d position, backPosition, velocity, cellCenter;
    Eigen::Array3i index;
    std::tie(i, j, k) = fluidFaces[it];

    if (_material[i][j][k] != Material::FluidMaterial::FLUID) {
      utemp[i][j][k] = _u[_currBuffer][i][j][k];
      vtemp[i][j][k] = _v[_currBuffer][i][j][k];
      wtemp[i][j][k] = _w[_currBuffer][i][j][k];
      continue;
    }
    double newVelocity = _u[_currBuffer][i][j][k];

    position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
               Eigen::Array3d(0, h[1] / 2, h[2] / 2);
    velocity[0] = _u[_currBuffer][i][j][k];
    velocity[1] = _interpolateVelocityV(position);
    velocity[2] = _interpolateVelocityW(position);

    backPosition = position - velocity * _dt;
    index = Eigen::floor(backPosition.cwiseQuotient(h)).cast<int>();

    if (velocity.matrix().norm() > _tolerance && index[0] >= 0 &&
        index[1] >= 0 && index[2] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1] && index[2] < _resolution[2]) {
      newVelocity = _interpolateVelocityU(backPosition);
    }
    utemp[i][j][k] = (newVelocity);
  }
  // -------------------------------------------
  fluidFaces.clear();
  for (auto id : _fluidCells) {
    Eigen::Array3i ijk = idToijk(id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];
    fluidFaces.emplace_back(std::make_tuple(i, j, k));
    if (j < _resolution[1] - 1 &&
        _material[i][j + 1][k] == Material::FluidMaterial::AIR) {
      fluidFaces.emplace_back(std::make_tuple(i, j + 1, k));
    }
  }
#pragma omp parallel for
  for (int it = 0; it < fluidFaces.size(); it++) {
    int i, j, k;
    Eigen::Array3d position, backPosition, velocity, cellCenter;
    Eigen::Array3i index;
    std::tie(i, j, k) = fluidFaces[it];

    double newVelocity = _v[_currBuffer][i][j][k];
    position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
               Eigen::Array3d(h[0] / 2, 0, h[2] / 2);
    velocity[0] = _interpolateVelocityU(position);
    velocity[1] = _v[_currBuffer][i][j][k];
    velocity[2] = _interpolateVelocityW(position);

    backPosition = position - velocity * _dt;
    index = Eigen::floor(backPosition.cwiseQuotient(h)).cast<int>();

    if ((velocity.matrix().norm() > _tolerance) && index[0] >= 0 &&
        index[1] >= 0 && index[2] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1] && index[2] < _resolution[2]) {
      newVelocity = _interpolateVelocityV(backPosition);
    }
    vtemp[i][j][k] = (newVelocity);
  }
  // -----------------------------------------------
  fluidFaces.clear();
  for (auto id : _fluidCells) {
    Eigen::Array3i ijk = idToijk(id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];
    fluidFaces.emplace_back(std::make_tuple(i, j, k));
    if (k < _resolution[2] &&
        _material[i][j][k + 1] == Material::FluidMaterial::AIR) {
      fluidFaces.emplace_back(std::make_tuple(i, j, k + 1));
    }
  }
#pragma omp parallel for
  for (int it = 0; it < fluidFaces.size(); it++) {
    int i, j, k;
    Eigen::Array3d position, backPosition, velocity, cellCenter;
    Eigen::Array3i index;
    std::tie(i, j, k) = fluidFaces[it];

    double newVelocity = _w[_currBuffer][i][j][k];

    position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
               Eigen::Array3d(h[0] / 2, h[1] / 2, 0);
    velocity[0] = _interpolateVelocityU(position);
    velocity[1] = _interpolateVelocityV(position);
    velocity[2] = _w[_currBuffer][i][j][k];

    backPosition = position - velocity * _dt;
    index = Eigen::floor(backPosition.cwiseQuotient(h)).cast<int>();

    if ((velocity.matrix().norm() > _tolerance) && index[0] >= 0 &&
        index[1] >= 0 && index[2] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1] && index[2] < _resolution[2]) {
      newVelocity = _interpolateVelocityW(backPosition);
    }
    wtemp[i][j][k] = (newVelocity);
  }

#pragma omp parallel for
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        if (utemp[i][j][k] > -1e7)
          _u[_currBuffer][i][j][k] = utemp[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        if (vtemp[i][j][k] > -1e7)
          _v[_currBuffer][i][j][k] = vtemp[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        if (wtemp[i][j][k] > -1e7)
          _w[_currBuffer][i][j][k] = wtemp[i][j][k];
      }
}

void RegularGrid3::setVelocity() {
  for (int i = 0; i < _resolution.x() + 1; i++) {
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _u[_currBuffer][i][j][k] = 0; // 1.0;
      }
  }
  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _v[_currBuffer][i][j][k] = 0; //-1.0;
      }
  }
  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        _w[_currBuffer][i][j][k] = 0; //-1.0;
      }
  }
}

int RegularGrid3::cellCount() {
  return _resolution.x() * _resolution.y() * _resolution.z();
}

int RegularGrid3::ijkToId(int i, int j, int k) {
  return k * _resolution.x() * _resolution.y() + j * _resolution.x() + i;
}

Eigen::Array3i RegularGrid3::idToijk(int cellId) {
  Eigen::Array3i index;
  index[2] = cellId / (_resolution.x() * _resolution.y());
  index[1] = (cellId % (_resolution.x() * _resolution.y())) / _resolution.x();
  index[0] = (cellId % (_resolution.x() * _resolution.y())) % _resolution.x();
  return index;
}

// TODO: Transfer writeFAceVelocity to and FileReader class
void RegularGrid3::writeFaceVelocity(const char *filename) {
  std::ofstream file;
  file.open(filename, std::ofstream::out | std::ofstream::trunc);
  if (!file.is_open()) {
    std::cerr << "RegularGrid3::writeFaceVelocity: unable to open " << filename
              << std::endl;
    return;
  }

  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x() + 1; i++) {
        file << _u[_currBuffer][i][j][k] << " ";
      }
      file << std::endl;
    }
    file << std::endl;
  }

  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y() + 1; j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        file << _v[_currBuffer][i][j][k] << " ";
      }
      file << std::endl;
    }
  file << std::endl;

  for (int k = 0; k < _resolution.z() + 1; k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        file << _w[_currBuffer][i][j][k] << " ";
      }
      file << std::endl;
    }
    file << std::endl;
  }
  file.close();
}

void RegularGrid3::readFaceVelocity(const char *filename) {
  std::ifstream file;
  file.open(filename);
  if (!file.is_open()) {
    std::cerr << "RegularGrid3::readFaceVelocity: unable to open " << filename
              << std::endl;
    return;
  }

  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x() + 1; i++) {
        file >> _u[_currBuffer][i][j][k];
      }
    }
  }

  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y() + 1; j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        file >> _v[_currBuffer][i][j][k];
      }
    }

  for (int k = 0; k < _resolution.z() + 1; k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        file >> _w[_currBuffer][i][j][k];
      }
    }
  }
  file.close();
}

void RegularGrid3::boundaryVelocities() {
  // Z velocities
  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      _w[_currBuffer][i][j][0] = 0.0;
      _w[_currBuffer][i][j][_resolution.z()] = 0.0;
    }
  }
  for (int k = 0; k < _resolution.z(); k++) {
    for (int i = 0; i < _resolution.x(); i++) {
      _v[_currBuffer][i][0][k] = 0.0;
      _v[_currBuffer][i][_resolution.y()][k] = 0.0;
    }
  }
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      _u[_currBuffer][0][j][k] = 0.0;
      _u[_currBuffer][_resolution.x()][j][k] = 0.0;
    }
  }
}

void RegularGrid3::addGravity() {
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 1; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        // if (_material[i][j][k] == Material::FluidMaterial::FLUID) {
        _v[_currBuffer][i][j][k] -= 9.81 * _dt;
        // if (j < _resolution[1] - 1 &&
        // _material[i][j + 1][k] == Material::FluidMaterial::AIR)
        // _v[_currBuffer][i][j + 1][k] -= 9.81 * _dt;
        // }extrapolateVel
      }
    }
  }
}

void RegularGrid3::addExternalForce(Eigen::Vector3d force) {
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];
    _u[_currBuffer][i][j][k] += force[0] * _dt;
    _v[_currBuffer][i][j][k] += force[1] * _dt;
    _w[_currBuffer][i][j][k] += force[2] * _dt;
  }
}

void RegularGrid3::extrapolateVelocity() {
  std::queue<int> processingCells;
  // TODO: Change processed cells to check phi instead of construct distance
  // vec
  Matrix3<int> processedCells;
  processedCells.changeSize(_resolution, 1e8);

  // Find first wavefront of surface fluid cells
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];
    // Check if is a air surface cell
    if (_material[i][j][k] == Material::FluidMaterial::FLUID) {
      processedCells[i][j][k] = 0;
      // Look for air neighbors
      // Add air cells to processing queue
      if (i > 0 && _material[i - 1][j][k] == Material::FluidMaterial::AIR) {
        processedCells[i - 1][j][k] = 1;
        processingCells.push(ijkToId(i - 1, j, k));
      }
      if (i < _resolution.x() - 1 &&
          _material[i + 1][j][k] == Material::FluidMaterial::AIR) {
        processedCells[i + 1][j][k] = 1;
        processingCells.push(ijkToId(i + 1, j, k));
      }
      if (j > 0 && _material[i][j - 1][k] == Material::FluidMaterial::AIR) {
        processedCells[i][j - 1][k] = 1;
        processingCells.push(ijkToId(i, j - 1, k));
      }
      if (j < _resolution.y() - 1 &&
          _material[i][j + 1][k] == Material::FluidMaterial::AIR) {
        processedCells[i][j + 1][k] = 1;
        processingCells.push(ijkToId(i, j + 1, k));
      }
      if (k > 0 && _material[i][j][k - 1] == Material::FluidMaterial::AIR) {
        processedCells[i][j][k - 1] = 1;
        processingCells.push(ijkToId(i, j, k - 1));
      }
      if (k < _resolution.y() - 1 &&
          _material[i][j][k + 1] == Material::FluidMaterial::AIR) {
        processedCells[i][j][k + 1] = 1;
        processingCells.push(ijkToId(i, j, k + 1));
      }
    }
  }

  // For each air cell, find its nearest neighbors from surface and take their
  // velocities values. In case when more than one is found, take the average
  while (!processingCells.empty()) {
    int currCell = processingCells.front();
    processingCells.pop();
    Eigen::Array3i ijk = idToijk(currCell);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    if (_material[i][j][k] == Material::FluidMaterial::AIR &&
        _hasOppositeNeighborsWithMaterial(currCell,
                                          Material::FluidMaterial::FLUID))
      continue;

    // if (processedCells[i][j][k] > 20)
    // continue;

    // Find the least distance
    int leastDistace = 1e6;
    std::vector<Eigen::Array3i> neighborCells;
    if (i > 0 && processedCells[i - 1][j][k] <= leastDistace) {
      if (processedCells[i - 1][j][k] < leastDistace) {
        leastDistace = processedCells[i - 1][j][k];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i - 1, j, k);
    }
    if (i < _resolution[0] - 1 && processedCells[i + 1][j][k] <= leastDistace) {
      if (processedCells[i + 1][j][k] < leastDistace) {
        leastDistace = processedCells[i + 1][j][k];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i + 1, j, k);
    }
    if (j > 0 && processedCells[i][j - 1][k] <= leastDistace) {
      if (processedCells[i][j - 1][k] < leastDistace) {
        leastDistace = processedCells[i][j - 1][k];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i, j - 1, k);
    }
    if (j < _resolution[1] - 1 && processedCells[i][j + 1][k] <= leastDistace) {
      if (processedCells[i][j + 1][k] < leastDistace) {
        leastDistace = processedCells[i][j + 1][k];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i, j + 1, k);
    }
    if (k > 0 && processedCells[i][j][k - 1] <= leastDistace) {
      if (processedCells[i][j][k - 1] < leastDistace) {
        leastDistace = processedCells[i][j][k - 1];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i, j, k - 1);
    }
    if (k < _resolution[2] - 1 && processedCells[i][j][k + 1] <= leastDistace) {
      if (processedCells[i][j][k + 1] < leastDistace) {
        leastDistace = processedCells[i][j][k + 1];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i, j, k + 1);
    }
    processedCells[i][j][k] = leastDistace + 1;
    // Find which faces need update
    // TODO: Change this to an enumerate
    bool faceNeedUpdate[] = {1, 1, 1,
                             1, 1, 1}; // LEFT RIGHT BOTTOM TOP BACK FRONT
    for (auto cell : neighborCells) {
      if (cell[0] < i)
        faceNeedUpdate[0] = false;
      else if (cell[0] > i)
        faceNeedUpdate[1] = false;
      if (cell[1] < j)
        faceNeedUpdate[2] = false;
      else if (cell[1] > j)
        faceNeedUpdate[3] = false;
      if (cell[2] < k)
        faceNeedUpdate[4] = false;
      else if (cell[2] > k)
        faceNeedUpdate[5] = false;
    }
    // Update the faces
    int nCells = neighborCells.size();
    if (faceNeedUpdate[0]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _u[_currBuffer][cell[0]][cell[1]][cell[2]];
      }
      _u[_currBuffer][i][j][k] = newVelocity / nCells;
      if (i > 0 && processedCells[i - 1][j][k] > 1e6) {
        processingCells.push(ijkToId(i - 1, j, k));
        processedCells[i - 1][j][k] = processedCells[i][j][k] + 1;
      }
    }
    if (faceNeedUpdate[1]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _u[_currBuffer][cell[0] + 1][cell[1]][cell[2]];
      }
      _u[_currBuffer][i + 1][j][k] = newVelocity / nCells;
      if (i < _resolution[0] - 1 && processedCells[i + 1][j][k] > 1e6) {
        processingCells.push(ijkToId(i + 1, j, k));
        processedCells[i + 1][j][k] = processedCells[i][j][k] + 1;
      }
    }
    if (faceNeedUpdate[2]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _v[_currBuffer][cell[0]][cell[1]][cell[2]];
      }
      _v[_currBuffer][i][j][k] = newVelocity / nCells;
      if (j > 0 && processedCells[i][j - 1][k] > 1e6) {
        processingCells.push(ijkToId(i, j - 1, k));
        processedCells[i][j - 1][k] = processedCells[i][j][k] + 1;
      }
    }
    if (faceNeedUpdate[3]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _v[_currBuffer][cell[0]][cell[1] + 1][cell[2]];
      }
      _v[_currBuffer][i][j + 1][k] = newVelocity / nCells;
      if (j < _resolution[1] - 1 && processedCells[i][j + 1][k] > 1e6) {
        processingCells.push(ijkToId(i, j + 1, k));
        processedCells[i][j + 1][k] = processedCells[i][j][k] + 1;
      }
    }
    if (faceNeedUpdate[4]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _w[_currBuffer][cell[0]][cell[1]][cell[2]];
      }
      _w[_currBuffer][i][j][k] = newVelocity / nCells;
      if (k > 0 && processedCells[i][j][k - 1] > 1e6) {
        processingCells.push(ijkToId(i, j, k - 1));
        processedCells[i][j][k - 1] = processedCells[i][j][k] + 1;
      }
    }
    if (faceNeedUpdate[5]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _w[_currBuffer][cell[0]][cell[1]][cell[2] + 1];
      }
      _w[_currBuffer][i][j][k + 1] = newVelocity / nCells;
      if (k < _resolution[2] - 1 && processedCells[i][j][k + 1] > 1e6) {
        processingCells.push(ijkToId(i, j, k + 1));
        processedCells[i][j][k + 1] = processedCells[i][j][k] + 1;
      }
    }
  }
}

void RegularGrid3::solvePressure() {

  // Compute velocity divergent over cell center
  // std::vector<double> divergent(cellCount(), 0.0);
  Eigen::VectorXd divergent, pressure;
  std::vector<Eigen::Triplet<double>> triplets;

  divergent = Eigen::VectorXd::Zero(cellCount());

  Timer timer;

  std::map<int, int> idMap;
  std::function<int(int)> _getMapId = [&](int cellId) {
    int index;
#pragma omp critical
    {
      auto it = idMap.find(cellId);
      if (it == idMap.end()) {
        int nextId = idMap.size();
        idMap[cellId] = nextId;
        index = nextId;
      } else
        index = it->second;
    }
    return index;
  };
// Solve pressure Poisson equation
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    if (_material[i][j][k] != Material::FluidMaterial::FLUID)
      continue;
    int id = _getMapId(_id);

    int validCells = 0;
    std::vector<Eigen::Triplet<double>> threadTriplet;

    if (i > 0) {
      validCells++;
      if (_material[i - 1][j][k] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i - 1, j, k)),
                                   -1 / (_h.x() * _h.x()));
    }
    if (i < _resolution.x() - 1) {
      validCells++;
      if (_material[i + 1][j][k] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i + 1, j, k)),
                                   -1 / (_h.x() * _h.x()));
    }
    if (j > 0) {
      validCells++;
      if (_material[i][j - 1][k] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j - 1, k)),
                                   -1 / (_h.x() * _h.x()));
    }
    if (j < _resolution.y() - 1) {
      validCells++;
      if (_material[i][j + 1][k] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j + 1, k)),
                                   -1 / (_h.x() * _h.x()));
    }
    if (k > 0) {
      validCells++;
      if (_material[i][j][k - 1] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j, k - 1)),
                                   -1 / (_h.x() * _h.x()));
    }
    if (k < _resolution.z() - 1) {
      validCells++;
      if (_material[i][j][k + 1] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j, k + 1)),
                                   -1 / (_h.x() * _h.x()));
    }

    threadTriplet.emplace_back(id, id, validCells / (_h.x() * _h.x()));
#pragma omp critical
    {
      triplets.insert(triplets.end(), threadTriplet.begin(),
                      threadTriplet.end());
    }

    divergent[id] = 0;
    divergent[id] -=
        (_u[_currBuffer][i + 1][j][k] - _u[_currBuffer][i][j][k]) / _h.x();
    divergent[id] -=
        (_v[_currBuffer][i][j + 1][k] - _v[_currBuffer][i][j][k]) / _h.y();
    divergent[id] -=
        (_w[_currBuffer][i][j][k + 1] - _w[_currBuffer][i][j][k]) / _h.z();
    divergent[id] /= _dt;
  }

  int nCells = idMap.size() + 1;
  triplets.emplace_back(nCells, nCells, 1);
  pressure = Eigen::VectorXd::Zero(nCells);
  divergent.conservativeResize(nCells);
  divergent[nCells] = 0.0;

  Eigen::SparseMatrix<double> pressureMatrix(nCells, nCells);
  pressureMatrix.setFromTriplets(triplets.begin(), triplets.end());
  timer.registerTime("Assembly");

  // SOlve pressure Poisson system
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(pressureMatrix);
  timer.registerTime("Prepare matrix");
  pressure = solver.solve(divergent);
  timer.registerTime("System solve");

// Correct velocity through pressure gradient
#pragma omp parallel for
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 1; i < _resolution.x(); i++) {
        double diffs =
            (pressure[ijkToId(i, j, k)] - pressure[ijkToId(i - 1, j, k)]) /
            _h.x();
        _u[_currBuffer][i][j][k] -= _dt * diffs;
      }
    }
  }
#pragma omp parallel for
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 1; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        double diffs =
            (pressure[ijkToId(i, j, k)] - pressure[ijkToId(i, j - 1, k)]) /
            _h.y();
        _v[_currBuffer][i][j][k] -= _dt * diffs;
      }
    }
  }
#pragma omp parallel for
  for (int k = 1; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        double diffs =
            (pressure[ijkToId(i, j, k)] - pressure[ijkToId(i, j, k - 1)]) /
            _h.z();
        _w[_currBuffer][i][j][k] -= _dt * diffs;
      }
    }
  }
  timer.registerTime("Gradient");
  timer.evaluateComponentsTime();
}

double RegularGrid3::_interpolateVelocityU(Eigen::Array3d position) {
  double _min, _max;
  _interpolateVelocityU(position, _min, _max);
}

double RegularGrid3::_interpolateVelocityU(Eigen::Array3d position,
                                           double &_min, double &_max) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution.x(), _resolution.y(), _resolution.z());
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  // Treating values that are outside domain

  if (position[1] < cellCenter[1])
    index[1]--;
  if (position[2] < cellCenter[2])
    index[2]--;

  std::vector<int> iCandidates, jCandidates, kCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
    kCandidates.push_back(index[2] + i);
  }
  std::vector<Eigen::Array3d> points;
  std::vector<double> values;

  for (auto w : kCandidates)
    for (auto v : jCandidates)
      for (auto u : iCandidates) {
        int x, y, z;
        points.emplace_back(u * h[0], (v + 0.5) * h[1], (w + 0.5) * h[2]);
        x = std::max(0, std::min(resolution[0], u));
        y = std::max(0, std::min(resolution[1] - 1, v));
        z = std::max(0, std::min(resolution[2] - 1, w));
        values.emplace_back(_u[_currBuffer][x][y][z]);
        _min = std::min(values.back(), _min);
        _max = std::max(values.back(), _max);
      }
  double velocity = Interpolator::tricubic(position, points, values);
  return velocity;
}

double RegularGrid3::_interpolateVelocityV(Eigen::Array3d position) {
  double _min, _max;
  _interpolateVelocityV(position, _min, _max);
}

double RegularGrid3::_interpolateVelocityV(Eigen::Array3d position,
                                           double &_min, double &_max) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution.x(), _resolution.y(), _resolution.z());
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  if (position[0] < cellCenter[0])
    index[0]--;
  if (position[2] < cellCenter[2])
    index[2]--;

  std::vector<int> iCandidates, jCandidates, kCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
    kCandidates.push_back(index[2] + i);
  }
  std::vector<Eigen::Array3d> points;
  std::vector<double> values;

  for (auto w : kCandidates)
    for (auto v : jCandidates)
      for (auto u : iCandidates) {
        int x, y, z;
        points.emplace_back((u + 0.5) * h[0], v * h[1], (w + 0.5) * h[2]);
        x = std::max(0, std::min(resolution[0] - 1, u));
        y = std::max(0, std::min(resolution[1], v));
        z = std::max(0, std::min(resolution[2] - 1, w));
        values.emplace_back(_v[_currBuffer][x][y][z]);
        _min = std::min(values.back(), _min);
        _max = std::max(values.back(), _max);
      }
  double velocity = Interpolator::tricubic(position, points, values);
  return velocity;
}

double RegularGrid3::_interpolateVelocityW(Eigen::Array3d position) {
  double _min, _max;
  _interpolateVelocityW(position, _min, _max);
}

double RegularGrid3::_interpolateVelocityW(Eigen::Array3d position,
                                           double &_min, double &_max) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution.x(), _resolution.y(), _resolution.z());
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  if (position[0] < cellCenter[0])
    index[0]--;
  if (position[1] < cellCenter[1])
    index[1]--;

  std::vector<int> iCandidates, jCandidates, kCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
    kCandidates.push_back(index[2] + i);
  }
  std::vector<Eigen::Array3d> points;
  std::vector<double> values;

  for (auto w : kCandidates)
    for (auto v : jCandidates)
      for (auto u : iCandidates) {
        int x, y, z;
        points.emplace_back((u + 0.5) * h[0], (v + 0.5) * h[1], w * h[2]);
        x = std::max(0, std::min(resolution[0] - 1, u));
        y = std::max(0, std::min(resolution[1] - 1, v));
        z = std::max(0, std::min(resolution[2], w));
        values.emplace_back(_w[_currBuffer][x][y][z]);
        _min = std::min(values.back(), _min);
        _max = std::max(values.back(), _max);
      }
  double velocity = Interpolator::tricubic(position, points, values);
  return velocity;
}

std::vector<Eigen::Array3i> RegularGrid3::cellFaces(Eigen::Array3i cell,
                                                    int coordinate) {
  std::vector<Eigen::Array3i> faces;

  faces.emplace_back(cell);
  if (coordinate == 0) {
    faces.emplace_back(cell + Eigen::Array3i(1, 0, 0));
  } else if (coordinate == 1) {
    faces.emplace_back(cell + Eigen::Array3i(0, 1, 0));
  } else if (coordinate == 2) {
    faces.emplace_back(cell + Eigen::Array3i(0, 0, 1));
  }

  return faces;
}

bool RegularGrid3::_hasOppositeNeighborsWithMaterial(
    int cellId, Material::FluidMaterial material) {
  Eigen::Array3i ijk = idToijk(cellId);
  int i, j, k;
  i = ijk[0];
  j = ijk[1];
  k = ijk[2];
  if (i > 0 && _material[i - 1][j][k] == material && i < _resolution[0] - 2 &&
      _material[i + 1][j][k] == material)
    return true;
  if (j > 0 && _material[i][j - 1][k] == material && j < _resolution[1] - 2 &&
      _material[i][j + 1][k] == material)
    return true;
  if (k > 0 && _material[i][j][k - 1] == material && k < _resolution[2] - 2 &&
      _material[i][j][k + 1] == material)
    return true;

  return false;
}

} // namespace Ramuh
