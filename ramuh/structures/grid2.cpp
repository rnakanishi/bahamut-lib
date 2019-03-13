#include <structures/grid2.h>
#include <utils/macros.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <queue>
#include <set>
#include <omp.h>

namespace Ramuh {

RegularGrid2::RegularGrid2() {
  setResolution(Vector2i(32, 32));
  _h = Vector2d(1. / 32);
  _domainSize = Vector2d(1.0);
  _dt = 0.01;
}

Vector2d RegularGrid2::domainSize() { return _domainSize; }

Vector2i RegularGrid2::resolution() { return _resolution; }

Vector2d RegularGrid2::h() { return _h; }

Vector2d RegularGrid2::gridSize() { return _domainSize; }

void RegularGrid2::setSize(Vector2d newSize) {
  _domainSize = newSize;
  setH(_domainSize / _resolution);
}

void RegularGrid2::setResolution(Vector2i newResolution) {
  _resolution = newResolution;
  setH(_domainSize / _resolution);

  _u.resize(_resolution.x() + 1);
  for (auto &row : _u) {
    row.resize(_resolution.y());
  }

  _v.resize(_resolution.x());
  for (auto &row : _v) {
    row.resize(_resolution.y() + 1);
  }

  _material.resize(_resolution.x());
  for (auto &row : _material) {
    row.resize(_resolution.y());
  }
}

void RegularGrid2::setH(Vector2d newH) { _h = newH; }

void RegularGrid2::advectGridVelocity() {
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
        vVel += _v[i - 1][j].y() + _v[i - 1][j + 1].y();
      }
      if (i < _resolution.x()) {
        vCount += 2;
        vVel += _v[i][j].y() + _v[i][j + 1].y();
      }
      utemp[i][j].x(_u[i][j].x());
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
        uVel += _u[i][j - 1].x() + _u[i + 1][j - 1].x();
      }
      if (j < _resolution.y()) {
        uCount += 2;
        uVel += _u[i][j].x() + _u[i + 1][j].x();
      }
      vtemp[i][j].x(uVel / uCount);
      vtemp[i][j].y(_v[i][j].y());
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
            _material[index.x()][index.y()];
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
        _u[i][j].x(utemp[i][j].x());
      else
        _u[i][j].x(newVel / distanceCount);
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
            _material[index.x()][index.y()];
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
        _v[i][j].y(vtemp[i][j].y());
      else
        _v[i][j].y(newVel / distanceCount);
    }

  // TODO: add convective term for z axis
}

void RegularGrid2::setVelocity() {

  for (int i = 0; i < _resolution.x() + 1; i++) {
    for (int j = 0; j < _resolution.y(); j++) {
      _u[i][j] = Vector2d(0, 0);
    }
  }

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y() + 1; j++) {
      _v[i][j] = Vector2d(0, 0);
    }
  }
}

int RegularGrid2::cellCount() { return _resolution.x() * _resolution.y(); }

int RegularGrid2::ijToId(int i, int j) { return j * _resolution.x() + i; }

Vector2i RegularGrid2::idToij(int id) { NOT_IMPLEMENTED(); }

void RegularGrid2::printFaceVelocity() {

  std::cerr << "==== u: \n";
  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _u[i][j] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << "==== v: \n";
  for (int j = 0; j < _resolution.y() + 1; j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      std::cerr << _v[i][j] << " ";
    }
    std::cerr << std::endl;
  }
}

void RegularGrid2::boundaryVelocities() {
  // Z velocities

  for (int i = 0; i < _resolution.x(); i++) {
    _v[i][0].y(0.0);
    _v[i][_resolution.y()].y(0.0);
  }
  for (int j = 0; j < _resolution.y(); j++) {
    _u[0][j].x(0.0);
    _u[_resolution.x()][j].x(0.0);
  }
}

void RegularGrid2::addGravity() {
  for (int j = 1; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      if (_material[i][j - 1] == Material::FluidMaterial::FLUID) {
        double vel = _v[i][j].y() - 9.81 * _dt;
        _v[i][j].y(vel);
      }
    }
  }
}

void RegularGrid2::extrapolateVelocity() {
  std::queue<Vector2i> processingCells;
  std::vector<std::vector<int>> processedCells;

  processedCells.resize(_resolution.x());
  for (auto &row : processedCells)
    row.resize(_resolution.y(), _resolution.x() * _resolution.y());

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.x(); j++) {
      // Extrapolate velocity from fluid to air cells
      // Check if is a air surface cell
      if (_material[i][j] == Material::FluidMaterial::FLUID) {
        processedCells[i][j] = 0;
        // Look for air neighbors
        // Add air cells to processing queue
        if (i > 0 && _material[i - 1][j] == Material::FluidMaterial::AIR) {
          processedCells[i - 1][j] = 1;
          processingCells.push(Vector2i(i - 1, j));
          _u[i - 1][j].x(_u[i][j].x());
          _v[i - 1][j].y(_v[i][j].y());
          _v[i - 1][j + 1].y(_v[i][j + 1].y());
        }
        if (i < _resolution.x() - 1 &&
            _material[i + 1][j] == Material::FluidMaterial::AIR) {
          processedCells[i + 1][j] = 1;
          processingCells.push(Vector2i(i + 1, j));
          _u[i + 2][j].x(_u[i + 1][j].x());
          _v[i + 1][j].y(_v[i][j].y());
          _v[i + 1][j + 1].y(_v[i][j + 1].y());
        }
        if (j > 0 && _material[i][j - 1] == Material::FluidMaterial::AIR) {
          processedCells[i][j - 1] = 1;
          processingCells.push(Vector2i(i, j - 1));
          _v[i][j - 1].y(_v[i][j].y());
          _u[i][j - 1].x(_u[i][j].x());
          _u[i + 1][j - 1].x(_u[i + 1][j].x());
        }
        if (j < _resolution.y() - 1 &&
            _material[i][j + 1] == Material::FluidMaterial::AIR) {
          processedCells[i][j + 1] = 1;
          processingCells.push(Vector2i(i, j + 1));
          _v[i][j + 2].y(_v[i][j + 1].y());
          _u[i][j + 1].x(_u[i][j].x());
          _u[i + 1][j + 1].x(_u[i + 1][j].x());
        }
      }
    }
  }
  while (!processingCells.empty()) {
    Vector2i currCell = processingCells.front();
    processingCells.pop();
    int i = currCell.x(), j = currCell.y();

    // Check for neighbor airs
    if (i > 0 && _material[i - 1][j] == Material::FluidMaterial::AIR) {
      if (processedCells[i - 1][j] > processedCells[i][j] + 1) {
        processedCells[i - 1][j] = processedCells[i][j] + 1;
        processingCells.push(Vector2i(i - 1, j));
        _u[i - 1][j].x(_u[i][j].x());
        _v[i - 1][j].y(_v[i][j].y());
        _v[i - 1][j + 1].y(_v[i][j + 1].y());
      }
    }
    if (i < _resolution.x() - 1 &&
        _material[i + 1][j] == Material::FluidMaterial::AIR) {
      if (processedCells[i + 1][j] > processedCells[i][j] + 1) {
        processedCells[i + 1][j] = processedCells[i][j] + 1;
        processingCells.push(Vector2i(i + 1, j));
        _u[i + 2][j].x(_u[i + 1][j].x());
        _v[i + 1][j].y(_v[i][j].y());
        _v[i + 1][j + 1].y(_v[i][j + 1].y());
      }
    }
    if (j > 0 && _material[i][j - 1] == Material::FluidMaterial::AIR) {
      if (processedCells[i][j - 1] > processedCells[i][j] + 1) {
        processedCells[i][j - 1] = processedCells[i][j] + 1;
        processingCells.push(Vector2i(i, j - 1));
        _v[i][j - 1].y(_v[i][j].y());
        _u[i][j - 1].x(_u[i][j].x());
        _u[i + 1][j - 1].x(_u[i + 1][j].x());
      }
    }
    if (j < _resolution.y() - 1 &&
        _material[i][j + 1] == Material::FluidMaterial::AIR) {
      if (processedCells[i][j + 1] > processedCells[i][j] + 1) {
        processedCells[i][j + 1] = processedCells[i][j] + 1;
        processingCells.push(Vector2i(i, j + 1));
        _v[i][j + 2].y(_v[i][j + 1].y());
        _u[i][j + 1].x(_u[i][j].x());
        _u[i + 1][j + 1].x(_u[i + 1][j].x());
      }
    }
  }
}

void RegularGrid2::solvePressure() {

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
  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      if (_material[i][j] != Material::FluidMaterial::FLUID)
        continue;

      int validCells = 0;
      int id = ijToId(i, j);
      std::vector<Eigen::Triplet<double>> threadTriplet;

      // TODO: impose second order boundary condition for pressure
      if (i > 0 || i < _resolution.x() - 1) {
        if (i > 0) {
          validCells++;
          if (_material[i - 1][j] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijToId(i - 1, j),
                                       -1 / (_h.x() * _h.x()));
        }
        if (i < _resolution.x() - 1) {
          validCells++;
          if (_material[i + 1][j] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijToId(i + 1, j),
                                       -1 / (_h.x() * _h.x()));
        }
      }
      if (j > 0 || j < _resolution.y() - 1) {
        if (j > 0) {
          validCells++;
          if (_material[i][j - 1] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijToId(i, j - 1),
                                       -1 / (_h.x() * _h.x()));
        }
        if (j < _resolution.y() - 1) {
          validCells++;
          if (_material[i][j + 1] != Material::FluidMaterial::AIR)
            threadTriplet.emplace_back(id, ijToId(i, j + 1),
                                       -1 / (_h.x() * _h.x()));
        }
      }
      threadTriplet.emplace_back(id, ijToId(i, j),
                                 validCells / (_h.x() * _h.x()));
#pragma omp critical
      {
        triplets.insert(triplets.end(), threadTriplet.begin(),
                        threadTriplet.end());
      }

      divergent[ijToId(i, j)] = 0;
      divergent[ijToId(i, j)] -= (_u[i + 1][j].x() - _u[i][j].x()) / _h.x();
      divergent[ijToId(i, j)] -= (_v[i][j + 1].y() - _v[i][j].y()) / _h.y();
      divergent[id] /= _dt;
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

  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 1; i < _resolution.x(); i++) {
      _u[i][j].x(_u[i][j].x() -
                 _dt * (pressure[ijToId(i, j)] - pressure[ijToId(i - 1, j)]) /
                     _h.x());
    }
  }

  for (int j = 1; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      _v[i][j].y(_v[i][j].y() -
                 _dt * (pressure[ijToId(i, j)] - pressure[ijToId(i, j - 1)]) /
                     _h.y());
    }
  }
}

} // namespace Ramuh
