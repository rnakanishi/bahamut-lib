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
#include <utils/timer.hpp>

namespace Ramuh {

RegularGrid3::RegularGrid3() {
  setResolution(Vector3i(32, 32, 32));
  _h = Vector3d(1. / 32);
  _domainSize = Vector3d(1.0);
  _dt = 0.01;
  _maxVelocity[0] = _maxVelocity[1] = _maxVelocity[2] = -1e8;
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

  double clamp[2]; // 0: min, 1: max
// Convective term for velocity U
#pragma omp parellel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
        Eigen::Vector3d velocity;
        Eigen::Array3i index;
        double newU = _u[i][j][k];
        // Find threshold values for clamping
        clamp[0] = 1e8;
        clamp[1] = -1e8;
        if (i > 0) {
          clamp[0] = std::min(clamp[0], _u[i - 1][j][k]);
          clamp[1] = std::max(clamp[1], _u[i - 1][j][k]);
        }
        if (i < _resolution[0]) {
          clamp[0] = std::min(clamp[0], _u[i + 1][j][k]);
          clamp[1] = std::max(clamp[1], _u[i + 1][j][k]);
        }
        if (j > 0) {
          clamp[0] = std::min(clamp[0], _u[i][j - 1][k]);
          clamp[1] = std::max(clamp[1], _u[i][j - 1][k]);
        }
        if (j < _resolution[1] - 1) {
          clamp[0] = std::min(clamp[0], _u[i][j + 1][k]);
          clamp[1] = std::max(clamp[1], _u[i][j + 1][k]);
        }
        if (k > 0) {
          clamp[0] = std::min(clamp[0], _u[i][j][k - 1]);
          clamp[1] = std::max(clamp[1], _u[i][j][k - 1]);
        }
        if (k < _resolution[2] - 1) {
          clamp[0] = std::min(clamp[0], _u[i][j][k + 1]);
          clamp[1] = std::max(clamp[1], _u[i][j][k + 1]);
        }

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(0, h[1] / 2, h[2] / 2);
        velocity[0] = _u[i][j][k];
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
          double error = 0.5 * (_u[i][j][k] - u_n1_hat);
          newU = u_n + error;
        } else {
          newU = _u[i][j][k];
        }
        utemp[i][j][k] = std::max(clamp[0], std::min(clamp[1], newU));
        //       }

        // // Convective term for velocity V
        // #pragma omp parellel for
        //   for (int i = 0; i < _resolution.x(); i++)
        //     for (int j = 1; j < _resolution.y(); j++)
        //       for (int k = 0; k < _resolution.z(); k++) {
        //         Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
        // Eigen::Vector3d velocity;
        // Eigen::Array3i index;
        double newV = _v[i][j][k];
        // Find threshold values for clamping
        clamp[0] = 1e8;
        clamp[1] = -1e8;
        if (i > 0) {
          clamp[0] = std::min(clamp[0], _v[i - 1][j][k]);
          clamp[1] = std::max(clamp[1], _v[i - 1][j][k]);
        }
        if (i < _resolution[0] - 1) {
          clamp[0] = std::min(clamp[0], _v[i + 1][j][k]);
          clamp[1] = std::max(clamp[1], _v[i + 1][j][k]);
        }
        if (j > 0) {
          clamp[0] = std::min(clamp[0], _v[i][j - 1][k]);
          clamp[1] = std::max(clamp[1], _v[i][j - 1][k]);
        }
        if (j < _resolution[1]) {
          clamp[0] = std::min(clamp[0], _v[i][j + 1][k]);
          clamp[1] = std::max(clamp[1], _v[i][j + 1][k]);
        }
        if (k > 0) {
          clamp[0] = std::min(clamp[0], _v[i][j][k - 1]);
          clamp[1] = std::max(clamp[1], _v[i][j][k - 1]);
        }
        if (k < _resolution[2] - 1) {
          clamp[0] = std::min(clamp[0], _v[i][j][k + 1]);
          clamp[1] = std::max(clamp[1], _v[i][j][k + 1]);
        }

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(h[0] / 2, 0, h[2] / 2);
        velocity[0] = _interpolateVelocityU(position);
        velocity[1] = _v[i][j][k];
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
          double error = 0.5 * (_v[i][j][k] - v_n1_hat);
          newV = v_n + error;
        } else {
          newV = _v[i][j][k];
        }
        vtemp[i][j][k] = std::max(clamp[0], std::min(clamp[1], newV));
        //       }

        // // Convective term for velocity W
        // #pragma omp parellel for
        //   for (int i = 0; i < _resolution.x(); i++)
        //     for (int j = 0; j < _resolution.y(); j++)
        //       for (int k = 1; k < _resolution.z(); k++) {
        //         Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
        // Eigen::Vector3d velocity;
        // Eigen::Array3i index;
        double newW = _w[i][j][k];
        // Find threshold values for clamping
        clamp[0] = 1e8;
        clamp[1] = -1e8;
        if (i > 0) {
          clamp[0] = std::min(clamp[0], _w[i - 1][j][k]);
          clamp[1] = std::max(clamp[1], _w[i - 1][j][k]);
        }
        if (i < _resolution[0] - 1) {
          clamp[0] = std::min(clamp[0], _w[i + 1][j][k]);
          clamp[1] = std::max(clamp[1], _w[i + 1][j][k]);
        }
        if (j > 0) {
          clamp[0] = std::min(clamp[0], _w[i][j - 1][k]);
          clamp[1] = std::max(clamp[1], _w[i][j - 1][k]);
        }
        if (j < _resolution[1] - 1) {
          clamp[0] = std::min(clamp[0], _w[i][j + 1][k]);
          clamp[1] = std::max(clamp[1], _w[i][j + 1][k]);
        }
        if (k > 0) {
          clamp[0] = std::min(clamp[0], _w[i][j][k - 1]);
          clamp[1] = std::max(clamp[1], _w[i][j][k - 1]);
        }
        if (k < _resolution[2]) {
          clamp[0] = std::min(clamp[0], _w[i][j][k + 1]);
          clamp[1] = std::max(clamp[1], _w[i][j][k + 1]);
        }

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(h[0] / 2, h[1] / 2, 0);
        velocity[0] = _interpolateVelocityU(position);
        velocity[1] = _interpolateVelocityV(position);
        velocity[2] = _w[i][j][k];
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
          double error = 0.5 * (_w[i][j][k] - w_n1_hat);
          newW = w_n + error;
        } else {
          newW = _w[i][j][k];
        }
        wtemp[i][j][k] = std::max(clamp[0], std::min(clamp[1], newW));
      }

#pragma omp barrier
#pragma omp parallel for
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _u[i][j][k] = utemp[i][j][k][0];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _v[i][j][k] = vtemp[i][j][k][1];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        _w[i][j][k] = wtemp[i][j][k][2];
      }
}

void RegularGrid3::advectGridVelocity() {
  // TODO: change to dynamic allocation
  Matrix3<double> utemp, vtemp, wtemp;
  Eigen::Array3d h(_h[0], _h[1], _h[2]);

  utemp.changeSize(
      Vector3i(_resolution.x() + 1, _resolution.y(), _resolution.z()));
  vtemp.changeSize(
      Vector3i(_resolution.x(), _resolution.y() + 1, _resolution.z()));
  wtemp.changeSize(
      Vector3i(_resolution.x(), _resolution.y(), _resolution.z() + 1));
#pragma omp parallel for
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        utemp[i][j][k] = _u[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        vtemp[i][j][k] = _v[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        wtemp[i][j][k] = _w[i][j][k];
      }

// Solve convection term (material derivative) - u velocities
#pragma omp parallel for
  for (int i = 0; i <= _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Eigen::Array3d position, backPosition, velocity, cellCenter;
        Eigen::Array3i index;
        double newVelocity = _u[i][j][k];

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(0, h[1] / 2, h[2] / 2);
        velocity[0] = _u[i][j][k];
        velocity[1] = _interpolateVelocityV(position);
        velocity[2] = _interpolateVelocityW(position);

        backPosition = position - velocity * _dt;
        index = Eigen::floor(backPosition.cwiseQuotient(h)).cast<int>();

        if (velocity.matrix().norm() > 1e-8 && index[0] >= 0 && index[1] >= 0 &&
            index[2] >= 0 && index[0] < _resolution[0] &&
            index[1] < _resolution[1] && index[2] < _resolution[2]) {
          newVelocity = _interpolateVelocityU(backPosition);
        }
        utemp[i][j][k] = (newVelocity);
      }

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j <= _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Eigen::Array3d position, backPosition, velocity, cellCenter;
        Eigen::Array3i index;
        double newVelocity = _v[i][j][k];

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(h[0] / 2, 0, h[2] / 2);
        velocity[0] = _interpolateVelocityU(position);
        velocity[1] = _v[i][j][k];
        velocity[2] = _interpolateVelocityW(position);

        backPosition = position - velocity * _dt;
        index = Eigen::floor(backPosition.cwiseQuotient(h)).cast<int>();

        if ((velocity.matrix().norm() > 1e-8) && index[0] >= 0 &&
            index[1] >= 0 && index[2] >= 0 && index[0] < _resolution[0] &&
            index[1] < _resolution[1] && index[2] < _resolution[2]) {
          newVelocity = _interpolateVelocityV(backPosition);
        }
        vtemp[i][j][k] = (newVelocity);
      }

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k <= _resolution.z(); k++) {
        Eigen::Array3d position, backPosition, velocity, cellCenter;
        Eigen::Array3i index;
        double newVelocity = _w[i][j][k];

        position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
                   Eigen::Array3d(h[0] / 2, h[1] / 2, 0);
        velocity[0] = _interpolateVelocityU(position);
        velocity[1] = _interpolateVelocityV(position);
        velocity[2] = _w[i][j][k];

        backPosition = position - velocity * _dt;
        index = Eigen::floor(backPosition.cwiseQuotient(h)).cast<int>();

        if ((velocity.matrix().norm() > 1e-8) && index[0] >= 0 &&
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
        _u[i][j][k] = utemp[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _v[i][j][k] = vtemp[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        _w[i][j][k] = wtemp[i][j][k];
      }
}

void RegularGrid3::setVelocity() {
  for (int i = 0; i < _resolution.x() + 1; i++) {
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _u[i][j][k] = 0; // 1.0;
      }
  }
  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _v[i][j][k] = 0; //-1.0;
      }
  }
  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        _w[i][j][k] = 0; //-1.0;
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
      std::cerr << _v[i][j][k] << " ";
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
      _w[i][j][0] = 0.0;
      _w[i][j][_resolution.z()] = 0.0;
    }
  }
  for (int k = 0; k < _resolution.z(); k++) {
    for (int i = 0; i < _resolution.x(); i++) {
      _v[i][0][k] = 0.0;
      _v[i][_resolution.y()][k] = 0.0;
    }
  }
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      _u[0][j][k] = 0.0;
      _u[_resolution.x()][j][k] = 0.0;
    }
  }
}

void RegularGrid3::addGravity() {
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 1; j <= _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        // if (_material[i][j][k] == Material::FluidMaterial::FLUID) {
        double vel = _v[i][j][k] - 9.81 * _dt;
        _v[i][j][k] = (vel);
        // }
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
    _u[i][j][k] += force[0] * _dt;
    _v[i][j][k] += force[1] * _dt;
    _w[i][j][k] += force[2] * _dt;
  }
}

void RegularGrid3::extrapolateVelocity() {
  std::queue<int> processingCells;
  // TODO: Change processed cells to check phi instead of construct distance vec
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
        // processedCells[i - 1][j][k] = 1;
        processingCells.push(ijkToId(i - 1, j, k));
      }
      if (i < _resolution.x() - 1 &&
          _material[i + 1][j][k] == Material::FluidMaterial::AIR) {
        // processedCells[i + 1][j][k] = 1;
        processingCells.push(ijkToId(i + 1, j, k));
      }
      if (j > 0 && _material[i][j - 1][k] == Material::FluidMaterial::AIR) {
        // processedCells[i][j - 1][k] = 1;
        processingCells.push(ijkToId(i, j - 1, k));
      }
      if (j < _resolution.y() - 1 &&
          _material[i][j + 1][k] == Material::FluidMaterial::AIR) {
        // processedCells[i][j + 1][k] = 1;
        processingCells.push(ijkToId(i, j + 1, k));
      }
      if (k > 0 && _material[i][j][k - 1] == Material::FluidMaterial::AIR) {
        // processedCells[i][j][k - 1] = 1;
        processingCells.push(ijkToId(i, j, k - 1));
      }
      if (k < _resolution.y() - 1 &&
          _material[i][j][k + 1] == Material::FluidMaterial::AIR) {
        // processedCells[i][j][k + 1] = 1;
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

    if (processedCells[i][j][k] < 1e7)
      continue;

    // Find the least distance
    int leastDistace = 1e8;
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
        newVelocity += _u[cell[0]][cell[1]][cell[2]];
      }
      _u[i][j][k] = newVelocity / nCells;
      if (i > 0 && processedCells[i - 1][j][k] > 1e6)
        processingCells.push(ijkToId(i - 1, j, k));
    }
    if (faceNeedUpdate[1]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _u[cell[0] + 1][cell[1]][cell[2]];
      }
      _u[i + 1][j][k] = newVelocity / nCells;
      if (i < _resolution[0] - 1 && processedCells[i + 1][j][k] > 1e6)
        processingCells.push(ijkToId(i + 1, j, k));
    }
    if (faceNeedUpdate[2]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _v[cell[0]][cell[1]][cell[2]];
      }
      _v[i][j][k] = newVelocity / nCells;
      if (j > 0 && processedCells[i][j - 1][k] > 1e6)
        processingCells.push(ijkToId(i, j - 1, k));
    }
    if (faceNeedUpdate[3]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _v[cell[0]][cell[1] + 1][cell[2]];
      }
      _v[i][j + 1][k] = newVelocity / nCells;
      if (j < _resolution[1] - 1 && processedCells[i][j + 1][k] > 1e6)
        processingCells.push(ijkToId(i, j + 1, k));
    }
    if (faceNeedUpdate[4]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _w[cell[0]][cell[1]][cell[2]];
      }
      _w[i][j][k] = newVelocity / nCells;
      if (k > 0 && processedCells[i][j][k - 1] > 1e6)
        processingCells.push(ijkToId(i, j, k - 1));
    }
    if (faceNeedUpdate[5]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _w[cell[0]][cell[1]][cell[2] + 1];
      }
      _w[i][j][k + 1] = newVelocity / nCells;
      if (k < _resolution[2] - 1 && processedCells[i][j][k + 1] > 1e6)
        processingCells.push(ijkToId(i, j, k + 1));
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

  Timer timer;

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
        divergent[id] -= (_u[i + 1][j][k] - _u[i][j][k]) / _h.x();
        divergent[id] -= (_v[i][j + 1][k] - _v[i][j][k]) / _h.y();
        divergent[id] -= (_w[i][j][k + 1] - _w[i][j][k]) / _h.z();
        divergent[id] /= _dt;
      }
    }
  }
  triplets.emplace_back(nCells, nCells, 1);

  pressureMatrix.setFromTriplets(triplets.begin(), triplets.end());
  timer.registerTime("Assembly");

  // SOlve pressure Poisson system
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;
  solver.compute(pressureMatrix);
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
        _u[i][j][k] -= _dt * diffs;
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
        _v[i][j][k] -= _dt * diffs;
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
        _w[i][j][k] -= _dt * diffs;
      }
    }
  }
  timer.registerTime("Gradient");
  timer.evaluateComponentsTime();
}

double RegularGrid3::_interpolateVelocityU(Eigen::Array3d position) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution.x(), _resolution.y(), _resolution.z());
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  std::vector<int> iCandidates, jCandidates, kCandidates;
  if (index[0] >= 0)
    iCandidates.push_back(index[0]);
  if (index[0] < _resolution[0] - 1)
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

  // Catmull-Rom like interpolation of u
  double clamp[2]; // 0: min, 1: max
  clamp[0] = 1e8;
  clamp[1] = -1e8;
  double distance = 0., distanceCount = 0.;
  double velocity = 0.0;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        if (u < 0 || u >= _resolution[0] || v < 0 || v >= _resolution[1] ||
            w < 0 || w >= _resolution[2])
          continue;
        Eigen::Array3d centerPosition =
            Eigen::Array3d(u, v, w) * h + Eigen::Array3d(0, h[1] / 2, h[2] / 2);
        distance = (position - centerPosition).matrix().norm();
        if (distance < 1e-6)
          return _u[u][v][w];
        distanceCount += 1. / distance;
        velocity += _u[u][v][w] / distance;
        clamp[0] = std::min(clamp[0], _u[u][v][w]);
        clamp[1] = std::max(clamp[1], _u[u][v][w]);
      }
  if (distanceCount == 0 || distanceCount > 1e8)
    throw(
        "RegularGrid3::interpolateVelocityU: distanceCount zero or infinity\n");
  return std::max(clamp[0], std::min(clamp[1], velocity / distanceCount));
  // return velocity / distanceCount;
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
  if (index[1] >= 0)
    jCandidates.push_back(index[1]);
  if (index[1] < _resolution[1] - 1)
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

  // Catmull-Rom like interpolation of v
  double clamp[2]; // 0: min, 1: max
  clamp[0] = 1e8;
  clamp[1] = -1e8;
  double distance = 0., distanceCount = 0.;
  double velocity = 0.0;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        if (u < 0 || u >= _resolution[0] || v < 0 || v >= _resolution[1] ||
            w < 0 || w >= _resolution[2])
          continue;
        Eigen::Array3d centerPosition =
            Eigen::Array3d(u, v, w) * h + Eigen::Array3d(h[0] / 2, 0, h[2] / 2);
        distance = (position - centerPosition).matrix().norm();
        if (distance < 1e-7)
          return _v[u][v][w];
        distanceCount += 1. / distance;
        velocity += _v[u][v][w] / distance;
        clamp[0] = std::min(clamp[0], _v[u][v][w]);
        clamp[1] = std::max(clamp[1], _v[u][v][w]);
      }
  if (distanceCount == 0 || distanceCount > 1e8)
    throw(
        "RegularGrid3::interpolateVelocityV: distanceCount zero or infinity\n");
  return std::max(clamp[0], std::min(clamp[1], velocity / distanceCount));
  // return velocity / distanceCount;
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
  if (index[2] >= 0)
    kCandidates.push_back(index[2]);
  if (index[2] < _resolution[2])
    kCandidates.push_back(index[2] + 1);
  if (position[0] > cellCenter[0] && index[0] < resolution[0] - 1)
    iCandidates.push_back(index[0] + 1);
  else if (position[0] < cellCenter[0] && index[0] > 0)
    iCandidates.push_back(index[0] - 1);
  if (position[1] > cellCenter[1] && index[1] < _resolution[1] - 1)
    jCandidates.push_back(index[1] + 1);
  else if (position[1] < cellCenter[1] && index[1] > 0)
    jCandidates.push_back(index[1] - 1);

  // Catmull-Rom like interpolation of w
  double clamp[2]; // 0: min, 1: max
  clamp[0] = 1e8;
  clamp[1] = -1e8;
  double distance = 0., distanceCount = 0.;
  double velocity = 0.0;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        if (u < 0 || u >= _resolution[0] || v < 0 || v >= _resolution[1] ||
            w < 0 || w >= _resolution[2])
          continue;
        Eigen::Array3d centerPosition =
            Eigen::Array3d(u, v, w) * h + Eigen::Array3d(h[0] / 2, h[1] / 2, 0);
        distance = (position - centerPosition).matrix().norm();
        if (distance < 1e-7)
          return _w[u][v][w];
        distanceCount += 1. / distance;
        velocity += _w[u][v][w] / distance;
        clamp[0] = std::min(clamp[0], _w[u][v][w]);
        clamp[1] = std::max(clamp[1], _w[u][v][w]);
      }
  if (distanceCount == 0 || distanceCount > 1e8)
    throw(
        "RegularGrid3::interpolateVelocityW: distanceCount zero or infinity\n");
  return std::max(clamp[0], std::min(clamp[1], velocity / distanceCount));
  // return velocity / distanceCount;
}

} // namespace Ramuh
