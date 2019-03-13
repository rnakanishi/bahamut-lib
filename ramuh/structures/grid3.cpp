#include <structures/grid3.h>
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

void RegularGrid3::advectGridVelocity() {
  // TODO: change to dynamic allocation
  std::vector<std::vector<Vector2d>> utemp;
  std::vector<std::vector<Vector2d>> vtemp;

  utemp.resize(_resolution.x() + 1);
  for (auto &row : utemp) {
    row.resize(_resolution.y());
  }

  vtemp.resize(_resolution.x());
  for (auto &row : vtemp) {
    row.resize(_resolution.y() + 1);
  }

// Interpolate other components
#pragma omp parallel for
  // x faces - all velocities components
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++) {
      double vVel = 0.0;
      int vCount = 0;

      if (i > 0) {
        vCount += 2;
        vVel += _v[i - 1][j][0].y() + _v[i - 1][j + 1][0].y();
      }
      if (i < _resolution.x()) {
        vCount += 2;
        vVel += _v[i][j][0].y() + _v[i][j + 1][0].y();
      }
      utemp[i][j].x(_u[i][j][0].x());
      utemp[i][j].y(vVel / vCount);
    }

#pragma omp parallel for
  // y faces - all velocities components
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y() + 1; j++) {
      double uVel = 0.0;
      int uCount = 0;

      if (j > 0) {
        uCount += 2;
        uVel += _u[i][j - 1][0].x() + _u[i + 1][j - 1][0].x();
      }
      if (j < _resolution.y()) {
        uCount += 2;
        uVel += _u[i][j][0].x() + _u[i + 1][j][0].x();
      }
      vtemp[i][j].x(uVel / uCount);
      vtemp[i][j].y(_v[i][j][0].y());
    }

// Solve convection term (material derivative) - u velocities
#pragma omp parallel for
  for (int i = 1; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++) {
      Vector2d h(_h.x(), _h.y());
      Vector2d backPosition,
          velPosition = Vector2d(i, j) * h + Vector2d(0, h.y() / 2);
      backPosition = velPosition - utemp[i][j] * _dt;

      // Velocity Linear interpolation
      Vector2i index;
      double newVel = 0.0;
      double distanceCount = 0.;
      double distance = 0.0;
      index.x(std::floor(backPosition.x() / h.x()));
      index.y(std::floor(backPosition.y() / h.y()));
      if (index >= Vector2i(0, 0) &&
          index < Vector2i(_resolution.x(), _resolution.y())) {
        Material::FluidMaterial centerMaterial =
            _material[index.x()][index.y()][0];
        Vector2d position;
        // double distance = (backPosition - position).length();
        // distanceCount += distance;

        if (index.y() >= 0) {
          if (index.x() >= 0) {
            // && _material[index.x()][index.y()][0] == centerMaterial) {
            position = index * h + (h * Vector2i(0, 0.5));
            distance = (backPosition - position).length();
            distanceCount += distance;
            newVel += distance * utemp[index.x()][index.y()].x();
          }
          if (index.x() < _resolution.x() - 1) {
            // && _material[index.x() + 1][index.y()][0] == centerMaterial) {
            position = (index + Vector2i(1, 0)) * h + (h * Vector2i(0, 0.5));
            distance = (backPosition - position).length();
            distanceCount += distance;
            newVel += distance * utemp[index.x() + 1][index.y()].x();
          }
        }
        if (index.y() < _resolution.y() - 1) {
          if (index.x() >= 0) {
            // && _material[index.x()][index.y() + 1][0] == centerMaterial) {
            position = (index + Vector2i(0, 1)) * h + (h * Vector2i(0, 0.5));
            distance = (backPosition - position).length();
            distanceCount += distance;
            newVel += distance * utemp[index.x()][index.y() + 1].x();
          }

          if (index.x() < _resolution.x() - 1) {
            // && _material[index.x() + 1][index.y() + 1][0] == centerMaterial)
            // {
            position = (index + Vector2i(1, 1)) * h + (h * Vector2i(0, 0.5));
            distance = (backPosition - position).length();
            distanceCount += distance;
            newVel += distance * utemp[index.x() + 1][index.y() + 1].x();
          }
        }
      } else {
        newVel = utemp[i][j].x();
        distanceCount = 1;
      }
      // TODO: Add verification for z axis
      if (distanceCount <= 1e-8)
        _u[i][j][0].x(utemp[i][j].x());
      else
        _u[i][j][0].x(newVel / distanceCount);
    }
// Solve convection term (material derivative) - v velocities
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 1; j < _resolution.y(); j++) {
      Vector2d h(_h.x(), _h.y());
      Vector2d backPosition,
          velPosition = Vector2d(i, j) * h + Vector2d(h.x() / 2, 0);
      backPosition = velPosition - vtemp[i][j] * _dt;

      // Velocity Linear interpolation
      Vector2i index;
      double newVel = 0.0;
      int distanceCount = 0;
      double distance = 0.0; //= (backPosition - position).length();
      index.x(std::floor(backPosition.x() / _h.x()));
      index.y(std::floor(backPosition.y() / _h.y()));
      // if (index.x() >= 0 || index.x() < _resolution.x() - 1) {
      if (index >= Vector2i(0, 0) &&
          index < Vector2i(_resolution.x(), _resolution.y())) {
        Material::FluidMaterial centerMaterial =
            _material[index.x()][index.y()][0];
        Vector2d position;
        // distanceCount += distance;

        if (index.x() >= 0) {
          if (index.y() >= 0) {
            // && _material[index.x()][index.y()][0] == centerMaterial) {
            position = index * h + (h * Vector2i(0.5, 0));
            distance = (backPosition - position).length();
            distanceCount += distance;
            newVel += distance * vtemp[index.x()][index.y()].y();
          }
          if (index.y() < _resolution.y() - 1) {
            // && _material[index.x()][index.y() + 1][0] == centerMaterial) {
            position = (index + Vector2i(0, 1)) * h + (h * Vector2i(0.5, 0));
            distance = (backPosition - position).length();
            distanceCount += distance;
            newVel += distance * vtemp[index.x()][index.y() + 1].y();
          }
        }
        if (index.x() < _resolution.x() - 1) {
          if (index.y() >= 0) {
            // && _material[index.x() + 1][index.y()][0] == centerMaterial) {
            position = (index + Vector2i(1, 0)) * h + (h * Vector2i(0.5, 0));
            distance = (backPosition - position).length();
            distanceCount += distance;
            newVel += distance * vtemp[index.x() + 1][index.y()].y();
          }
          if (index.y() < _resolution.y() - 1) {
            // && _material[index.x()][index.y() + 1][0] == centerMaterial) {
            position = (index + Vector2i(1, 1)) * h + (h * Vector2i(0.5, 0));
            distance = (backPosition - position).length();
            distanceCount += distance;
            newVel += distance * vtemp[index.x() + 1][index.y() + 1].y();
          }
        }
      } else {
        newVel = vtemp[i][j].y();
        distanceCount = 1;
      }
      if (distanceCount == 0)
        _v[i][j][0].y(vtemp[i][j].y());
      else
        _v[i][j][0].y(newVel / distanceCount);
    }

  // TODO: add convective term for z axis
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

  std::cerr << "==== u: \n";
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x() + 1; i++) {
        std::cerr << _u[i][j][0] << " ";
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }
  std::cerr << "==== v: \n";
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y() + 1; j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        std::cerr << _v[i][j][0] << " ";
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }
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
        if (_material[i][j - 1][k] == Material::FluidMaterial::FLUID) {
          double vel = _v[i][j][k].y() - 9.81 * _dt;
          _v[i][j][k].y(vel);
        }
      }
    }
  }
}

void RegularGrid3::extrapolateVelocity() {
  std::queue<Vector3i> processingCells;
  std::vector<std::vector<int>> processedCells;

  processedCells.resize(_resolution.x());
  for (auto &row : processedCells)
    row.resize(_resolution.y(), _resolution.x() * _resolution.y());

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.x(); j++) {
      // Extrapolate velocity from fluid to air cells
      // Check if is a air surface cell
      if (_material[i][j][0] == Material::FluidMaterial::FLUID) {
        processedCells[i][j] = 0;
        // Look for air neighbors
        // Add air cells to processing queue
        if (i > 0 && _material[i - 1][j][0] == Material::FluidMaterial::AIR) {
          processedCells[i - 1][j] = 1;
          processingCells.push(Vector3i(i - 1, j, 0));
          _u[i - 1][j][0].x(_u[i][j][0].x());
          _v[i - 1][j][0].y(_v[i][j][0].y());
          _v[i - 1][j + 1][0].y(_v[i][j + 1][0].y());
        }
        if (i < _resolution.x() - 1 &&
            _material[i + 1][j][0] == Material::FluidMaterial::AIR) {
          processedCells[i + 1][j] = 1;
          processingCells.push(Vector3i(i + 1, j, 0));
          _u[i + 2][j][0].x(_u[i + 1][j][0].x());
          _v[i + 1][j][0].y(_v[i][j][0].y());
          _v[i + 1][j + 1][0].y(_v[i][j + 1][0].y());
        }
        if (j > 0 && _material[i][j - 1][0] == Material::FluidMaterial::AIR) {
          processedCells[i][j - 1] = 1;
          processingCells.push(Vector3i(i, j - 1, 0));
          _v[i][j - 1][0].y(_v[i][j][0].y());
          _u[i][j - 1][0].x(_u[i][j][0].x());
          _u[i + 1][j - 1][0].x(_u[i + 1][j][0].x());
        }
        if (j < _resolution.y() - 1 &&
            _material[i][j + 1][0] == Material::FluidMaterial::AIR) {
          processedCells[i][j + 1] = 1;
          processingCells.push(Vector3i(i, j + 1, 0));
          _v[i][j + 2][0].y(_v[i][j + 1][0].y());
          _u[i][j + 1][0].x(_u[i][j][0].x());
          _u[i + 1][j + 1][0].x(_u[i + 1][j][0].x());
        }
      }
    }
  }
  while (!processingCells.empty()) {
    Vector3i currCell = processingCells.front();
    processingCells.pop();
    int i = currCell.x(), j = currCell.y();

    // Check for neighbor airs
    if (i > 0 && _material[i - 1][j][0] == Material::FluidMaterial::AIR) {
      if (processedCells[i - 1][j] > processedCells[i][j] + 1) {
        processedCells[i - 1][j] = processedCells[i][j] + 1;
        processingCells.push(Vector3i(i - 1, j, 0));
        _u[i - 1][j][0].x(_u[i][j][0].x());
        _v[i - 1][j][0].y(_v[i][j][0].y());
        _v[i - 1][j + 1][0].y(_v[i][j + 1][0].y());
      }
    }
    if (i < _resolution.x() - 1 &&
        _material[i + 1][j][0] == Material::FluidMaterial::AIR) {
      if (processedCells[i + 1][j] > processedCells[i][j] + 1) {
        processedCells[i + 1][j] = processedCells[i][j] + 1;
        processingCells.push(Vector3i(i + 1, j, 0));
        _u[i + 2][j][0].x(_u[i + 1][j][0].x());
        _v[i + 1][j][0].y(_v[i][j][0].y());
        _v[i + 1][j + 1][0].y(_v[i][j + 1][0].y());
      }
    }
    if (j > 0 && _material[i][j - 1][0] == Material::FluidMaterial::AIR) {
      if (processedCells[i][j - 1] > processedCells[i][j] + 1) {
        processedCells[i][j - 1] = processedCells[i][j] + 1;
        processingCells.push(Vector3i(i, j - 1, 0));
        _v[i][j - 1][0].y(_v[i][j][0].y());
        _u[i][j - 1][0].x(_u[i][j][0].x());
        _u[i + 1][j - 1][0].x(_u[i + 1][j][0].x());
      }
    }
    if (j < _resolution.y() - 1 &&
        _material[i][j + 1][0] == Material::FluidMaterial::AIR) {
      if (processedCells[i][j + 1] > processedCells[i][j] + 1) {
        processedCells[i][j + 1] = processedCells[i][j] + 1;
        processingCells.push(Vector3i(i, j + 1, 0));
        _v[i][j + 2][0].y(_v[i][j + 1][0].y());
        _u[i][j + 1][0].x(_u[i][j][0].x());
        _u[i + 1][j + 1][0].x(_u[i + 1][j][0].x());
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

#pragma omp parallel
  // Solve pressure Poisson equation
  {
    for (int k = 0; k < _resolution.z(); k++) {
#pragma omp for
      for (int j = 0; j < _resolution.y(); j++) {
        for (int i = 0; i < _resolution.x(); i++) {
          if (_material[i][j][k] != Material::FluidMaterial::FLUID)
            continue;

          int validCells = 0;
          int id = ijkToId(i, j, k);
          std::vector<Eigen::Triplet<double>> threadTriplet;

          if (i > 0 || i < _resolution.x() - 1) {
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
          }
          if (j > 0 || j < _resolution.y() - 1) {
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
          }
          if (k > 0 || k < _resolution.z() - 1) {
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
          }
          threadTriplet.emplace_back(id, ijkToId(i, j, k),
                                     validCells / (_h.x() * _h.x()));
#pragma omp critical
          {
            triplets.insert(triplets.end(), threadTriplet.begin(),
                            threadTriplet.end());
          }

          divergent[ijkToId(i, j, k)] = 0;
          divergent[ijkToId(i, j, k)] -=
              (_u[i + 1][j][k].x() - _u[i][j][k].x()) / _h.x();
          divergent[ijkToId(i, j, k)] -=
              (_v[i][j + 1][k].y() - _v[i][j][k].y()) / _h.y();
          divergent[ijkToId(i, j, k)] -=
              (_w[i][j][k + 1].z() - _w[i][j][k].z()) / _h.z();
          divergent[id] /= _dt;
        }
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
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 1; i < _resolution.x(); i++) {
        _u[i][j][k].x(_u[i][j][k].x() - _dt *
                                            (pressure[ijkToId(i, j, k)] -
                                             pressure[ijkToId(i - 1, j, k)]) /
                                            _h.x());
      }
    }
  }
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 1; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        _v[i][j][k].y(_v[i][j][k].y() - _dt *
                                            (pressure[ijkToId(i, j, k)] -
                                             pressure[ijkToId(i, j - 1, k)]) /
                                            _h.y());
      }
    }
  }
  for (int k = 1; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y(); j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        _w[i][j][k].z(_w[i][j][k].z() - _dt *
                                            (pressure[ijkToId(i, j, k)] -
                                             pressure[ijkToId(i, j, k - 1)]) /
                                            _h.z());
      }
    }
  }
}

} // namespace Ramuh
