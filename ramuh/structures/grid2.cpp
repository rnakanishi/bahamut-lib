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
      Vector2d position, backPosition, velocity, h(_h.x(), _h.y()), cellCenter;
      Vector2i index;
      Vector2d newVelocity = _u[i][j];
      double distanceCount = 0.0, distance = 0.0;

      velocity = utemp[i][j];
      position = Vector2d(i, j) * h + Vector2d(0, h.y() / 2);
      backPosition = position - velocity * _dt;
      index.x(std::floor(backPosition.x() / h.x()));
      index.y(std::floor(backPosition.y() / h.y()));

      if (velocity.length() > 1e-8 && index >= Vector2i(0, 0) &&
          index < Vector2i(_resolution.x(), _resolution.y())) {
        // Inside simulation domain
        cellCenter = Vector2d(index.x(), index.y()) * h + h / 2.0;
        std::vector<int> iCandidates, jCandidates;
        iCandidates.push_back(index.x());
        iCandidates.push_back(index.x() + 1);
        jCandidates.push_back(index.y());
        if (backPosition.y() > cellCenter.y() &&
            index.y() < _resolution.y() - 1)
          jCandidates.push_back(index.y() + 1);
        else if (backPosition.y() < cellCenter.y() && index.y() > 0)
          jCandidates.push_back(index.y() - 1);

        newVelocity = 0.0;
        for (auto u : iCandidates)
          for (auto v : jCandidates) {
            position = Vector2i(u, v) * h + Vector2d(0, h.y() / 2);
            distance = (backPosition - position).length();
            distanceCount += 1. / distance;
            newVelocity = newVelocity + (_u[u][v] * (1. / distance));
          }
        newVelocity = newVelocity / distanceCount;
      }
      _u[i][j] = newVelocity;
    }

// Solve convection term (material derivative) - v velocities
#pragma omp for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 1; j < _resolution.y(); j++) {
      Vector2d position, backPosition, velocity, h(_h.x(), _h.y()), cellCenter;
      Vector2i index;
      double newVelocity = _v[i][j].y();
      double distanceCount = 0.0, distance = 0.0;

      velocity = vtemp[i][j];
      position = Vector2d(i, j) * h + Vector2d(h.x() / 2, 0);
      backPosition = position - velocity * _dt;
      index.x(std::floor(backPosition.x() / h.x()));
      index.y(std::floor(backPosition.y() / h.y()));

      if (velocity.length() > 1e-8 && index >= Vector2i(0, 0) &&
          index < Vector2i(_resolution.x(), _resolution.y())) {
        // Inside simulation domain
        cellCenter = Vector2d(index.x(), index.y()) * h + h / 2.0;
        std::vector<int> iCandidates, jCandidates;
        iCandidates.push_back(index.x());
        jCandidates.push_back(index.y());
        jCandidates.push_back(index.y() + 1);

        if (backPosition.x() > cellCenter.x() &&
            index.x() < _resolution.x() - 1)
          iCandidates.push_back(index.x() + 1);
        else if (backPosition.x() < cellCenter.x() && index.x() > 0)
          iCandidates.push_back(index.x() - 1);

        newVelocity = 0.0;
        for (auto u : iCandidates)
          for (auto v : jCandidates) {
            position = Vector2i(u, v) * h + Vector2d(h.x() / 2, 0);
            distance = (backPosition - position).length();
            distanceCount += 1. / distance;
            newVelocity += (vtemp[u][v].y() * (1. / distance));
          }
        newVelocity /= distanceCount;
      }
      _v[i][j].y(newVelocity);
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
  for (int j = _resolution.y() - 1; j >= 0; j--) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _u[i][j].x() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << "==== v: \n";
  for (int j = _resolution.y() - 1; j >= 0; j--) {
    for (int i = 0; i < _resolution.x(); i++) {
      std::cerr << _v[i][j].y() << " ";
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
      // if (_material[i][j - 1] == Material::FluidMaterial::FLUID) {
      double vel = _v[i][j].y() - 9.81 * _dt;
      _v[i][j].y(vel);
      // }
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
  // Assemble Laplacian matrix
  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      if (_material[i][j] != Material::FluidMaterial::FLUID)
        continue;

      int validCells = 0;
      int id = ijToId(i, j);
      std::vector<Eigen::Triplet<double>> threadTriplet;

      // TODO: impose second order boundary condition for pressure
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
  // std::cerr << divergent.transpose() << std::endl;
  // std::cerr << pressure.transpose() << std::endl;
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
