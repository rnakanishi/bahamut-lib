#include <structures/grid2.h>
#include <utils/macros.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <queue>
#include <set>
#include <omp.h>
#include <utils/timer.hpp>
#include <blas/interpolator.h>

namespace Ramuh {

RegularGrid2::RegularGrid2() {
  setResolution(Eigen::Array2i(32, 32));
  _h = Eigen::Array2d(1. / 32, 1. / 32);
  _domainSize = Eigen::Array2d(1.0);

  _ellapsedDt = 0.0;
  _originalDt = _dt = 1.0 / 60;
  _tolerance = 1e-12;
}

Eigen::Array2d RegularGrid2::domainSize() { return _domainSize; }

Eigen::Array2i RegularGrid2::resolution() { return _resolution; }

Eigen::Array2d RegularGrid2::h() { return _h; }

Eigen::Array2d RegularGrid2::gridSize() { return _domainSize; }

void RegularGrid2::setDt(double dt) {
  _originalDt = dt;
  _dt = dt;
}

void RegularGrid2::setSize(Eigen::Array2d newSize) {
  _domainSize = newSize;
  setH(_domainSize.cwiseQuotient(_resolution.cast<double>()));
}

void RegularGrid2::setResolution(Eigen::Array2i newResolution) {
  _resolution = newResolution;
  setH(_domainSize.cwiseQuotient(_resolution.cast<double>()));

  _u.resize(_resolution[0] + 1);
  _uFaceMaterial.resize(_resolution[0] + 1);
  for (int i = 0; i < _resolution[0] + 1; i++) {
    _u[i].resize(_resolution[1]);
    _uFaceMaterial[i].resize(_resolution[1]);
  }

  _v.resize(_resolution[0]);
  _vFaceMaterial.resize(_resolution[0]);
  for (int i = 0; i < _resolution[0] + 1; i++) {
    _u[i].resize(_resolution[1]);
    _uFaceMaterial[i].resize(_resolution[1]);
  }

  _material.resize(_resolution[0]);
  for (auto &row : _material) {
    row.resize(_resolution[1]);
  }
}

void RegularGrid2::setH(Eigen::Array2d newH) { _h = newH; }

bool RegularGrid2::advanceTime() {
  _ellapsedDt += _dt;
  if (_ellapsedDt < _originalDt)
    return false;
  _ellapsedDt = 0.0;
  _dt = _originalDt;
  return true;
}

void RegularGrid2::cfl() {
  // Find biggest velocity

  Eigen::Vector2d maxVel = Eigen::Vector2d(0.0, 0.0);
  for (int id = 0; id < cellCount(); id++) {
    Eigen::Array2i ij = idToij(id);
    int i, j;
    i = ij[0];
    j = ij[1];

    Eigen::Vector2d vel;
    vel[0] = (_u[i][j] + _u[i + 1][j]) / 2.0;
    vel[1] = (_v[i][j] + _v[i][j + 1]) / 2.0;

    if (maxVel.norm() < vel.norm() * 1.05)
      maxVel = vel;
  }

  // Check if cfl condition applies
  // Half timestep if so
  if (maxVel.norm() * (_originalDt - _ellapsedDt) > 3 * _h[0]) {
    _dt = _dt / 2;
  }
}

void RegularGrid2::macCormackVelocityAdvection() {
  std::vector<std::vector<double>> utemp, vtemp;
  utemp.resize(_resolution.x() + 1);
  for (auto &row : utemp) {
    row.resize(_resolution[1]);
  }
  vtemp.resize(_resolution.x() + 1);
  for (auto &row : vtemp) {
    row.resize(_resolution[1]);
  }
  double clamp[2]; // 0: min, 1: max
// Convective term for velocity U
#pragma omp parellel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution[1]; j++) {
      Eigen::Array2d position, h(_h[0], _h[1]);
      Eigen::Vector2d velocity;
      Eigen::Array2i index;

      if (_uFaceMaterial[i][j] != Material::FluidMaterial::FLUID) {
        utemp[i][j] = _u[i][j];
      } else {
        double newU = _u[i][j];

        position =
            Eigen::Array2d(i, j).cwiseProduct(h) + Eigen::Array2d(0, h[1] / 2);
        velocity[0] = _u[i][j];
        velocity[1] = _interpolateVelocityV(position);
        position -= velocity.array() * _dt;
        index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();

        if ((velocity.norm() > _tolerance) && index[0] >= 0 && index[1] >= 0 &&
            index[0] < _resolution[0] && index[1] < _resolution[1]) {
          double u_n;
          u_n = _interpolateVelocityU(position, clamp[0], clamp[1]);
          velocity[0] = _interpolateVelocityU(position);
          velocity[1] = _interpolateVelocityV(position);
          position += velocity.array() * _dt;
          double u_n1_hat = _interpolateVelocityU(position);
          double error = 0.5 * (_u[i][j] - u_n1_hat);
          newU = u_n + error;
        } else {
          newU = _u[i][j];
        }
        utemp[i][j] = std::max(clamp[0], std::min(clamp[1], newU));
      }
      // utemp[i][j] = newU;
      //       }

      // // Convective term for velocity V
      if (_uFaceMaterial[i][j] != Material::FluidMaterial::FLUID) {
        vtemp[i][j] = _v[i][j];
      } else {
        double newV = _v[i][j];

        position =
            Eigen::Array2d(i, j).cwiseProduct(h) + Eigen::Array2d(h[0] / 2, 0);
        velocity[0] = _interpolateVelocityU(position);
        velocity[1] = _v[i][j];
        position -= velocity.array() * _dt;
        index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();

        if ((velocity.norm() > _tolerance) && index[0] >= 0 && index[1] >= 0 &&
            index[0] < _resolution[0] && index[1] < _resolution[1]) {
          double v_n;
          v_n = _interpolateVelocityV(position, clamp[0], clamp[1]);
          velocity[0] = _interpolateVelocityU(position);
          velocity[1] = _interpolateVelocityV(position);
          position += velocity.array() * _dt;
          double v_n1_hat = _interpolateVelocityV(position);
          double error = 0.5 * (_v[i][j] - v_n1_hat);
          newV = v_n + error;
        } else {
          newV = _v[i][j];
        }
        vtemp[i][j] = std::max(clamp[0], std::min(clamp[1], newV));
      }
      // vtemp[i][j] = newV;
      //       }
    }

#pragma omp barrier
#pragma omp parallel for
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution[1]; j++) {
      _u[i][j] = utemp[i][j];
    }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution[1] + 1; j++) {
      _v[i][j] = vtemp[i][j];
    }
}

void RegularGrid2::advectGridVelocity() {
  // TODO: change to dynamic allocation
  std::vector<std::vector<double>> utemp;
  std::vector<std::vector<double>> vtemp;

  utemp.resize(_resolution.x() + 1);
  for (auto &row : utemp) {
    row.resize(_resolution[1]);
  }

  vtemp.resize(_resolution.x());
  for (auto &row : vtemp) {
    row.resize(_resolution[1] + 1);
  }
  std::vector<std::tuple<int, int>> fluidFaces;
  for (auto id : _fluidCells) {
    Eigen::Array2i ij = idToij(id);
    int i, j;
    i = ij[0];
    j = ij[1];
    fluidFaces.emplace_back(std::make_tuple(i, j));
    if (i < _resolution[0] - 1 &&
        _material[i + 1][j] == Material::FluidMaterial::AIR) {
      fluidFaces.emplace_back(std::make_tuple(i + 1, j));
    }
  }
  // Interpolate other components
  Eigen::Array2i ij;
#pragma omp parallel for
  for (int it = 0; it < fluidFaces.size(); it++) {
    int i, j;
    Eigen::Array2d position, backPosition, velocity, cellCenter;
    Eigen::Array2i index;
    std::tie(i, j) = fluidFaces[it];

    if (_material[i][j] != Material::FluidMaterial::FLUID) {
      utemp[i][j] = _u[i][j];
      vtemp[i][j] = _v[i][j];
      continue;
    }
    double newVelocity = _u[i][j];

    position =
        Eigen::Array2d(i, j).cwiseProduct(_h) + Eigen::Array2d(0.0, _h[1] / 2);
    velocity[0] = _u[i][j];
    velocity[1] = _interpolateVelocityV(position);

    position = position - velocity * _dt;
    index = Eigen::floor(position.cwiseQuotient(_h)).cast<int>();

    if (velocity.matrix().norm() > _tolerance && index[0] >= 0 &&
        index[1] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1]) {
      newVelocity = _interpolateVelocityU(backPosition);
    }
    utemp[i][j] = (newVelocity);
  }
  //--------------------------------------------------------
  fluidFaces.clear();
  for (auto id : _fluidCells) {
    Eigen::Array2i ij = idToij(id);
    int i, j;
    i = ij[0];
    j = ij[1];
    fluidFaces.emplace_back(std::make_tuple(i, j));
    if (j < _resolution[1] - 1 &&
        _material[i][j + 1] == Material::FluidMaterial::AIR) {
      fluidFaces.emplace_back(std::make_tuple(i, j + 1));
    }
  }
#pragma omp parallel for
  for (int it = 0; it < fluidFaces.size(); it++) {
    int i, j;
    Eigen::Array2d position, backPosition, velocity, cellCenter;
    Eigen::Array2i index;
    std::tie(i, j) = fluidFaces[it];

    double newVelocity = _v[i][j];
    position =
        Eigen::Array2d(i, j).cwiseProduct(_h) + Eigen::Array2d(_h[0] / 2, 0.0);
    velocity[0] = _interpolateVelocityU(position);
    velocity[1] = _v[i][j];

    position = position - velocity * _dt;
    index = Eigen::floor(position.cwiseQuotient(_h)).cast<int>();

    if ((velocity.matrix().norm() > _tolerance) && index[0] >= 0 &&
        index[1] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1]) {
      newVelocity = _interpolateVelocityV(backPosition);
    }
    vtemp[i][j] = (newVelocity);
  }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution[1]; j++) {
      if (utemp[i][j] > -1e7)
        _u[i][j] = utemp[i][j];
    }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution[1] + 1; j++) {
      if (vtemp[i][j] > -1e7)
        _v[i][j] = vtemp[i][j];
    }
}

void RegularGrid2::setVelocity() {

  for (int i = 0; i < _resolution.x() + 1; i++) {
    for (int j = 0; j < _resolution[1]; j++) {
      _u[i][j] = 0.0;
    }
  }

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution[1] + 1; j++) {
      _v[i][j] = 0.0;
    }
  }
}

int RegularGrid2::cellCount() { return _resolution.x() * _resolution[1]; }

int RegularGrid2::ijToId(int i, int j) { return j * _resolution.x() + i; }

Eigen::Array2i RegularGrid2::idToij(int id) {
  Eigen::Array2i index;
  index[1] = id / (_resolution[0]);
  index[0] = id % _resolution[0];
  return index;
}

void RegularGrid2::printFaceVelocity() {

  std::cerr << "==== u: \n";
  for (int j = _resolution[1] - 1; j >= 0; j--) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _u[i][j] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << "==== v: \n";
  for (int j = _resolution[1] - 1; j >= 0; j--) {
    for (int i = 0; i < _resolution.x(); i++) {
      std::cerr << _v[i][j] << " ";
    }
    std::cerr << std::endl;
  }
}

void RegularGrid2::boundaryVelocities() {
  // Z velocities
  for (int i = 0; i < _resolution.x(); i++) {
    _v[i][0] = 0.0;
    _v[i][_resolution[1]] = 0.0;
  }
  for (int j = 0; j < _resolution[1]; j++) {
    _u[0][j] = 0.0;
    _u[_resolution.x()][j] = 0.0;
  }
}

void RegularGrid2::addGravity() {
  for (int j = 1; j < _resolution[1]; j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      // if (_material[i][j - 1] == Material::FluidMaterial::FLUID) {
      _v[i][j] = _v[i][j] - 9.81 * _dt;
      // }
    }
  }
}

void RegularGrid2::extrapolateVelocity() {
  std::queue<int> processingCells;
  // TODO: Change processed cells to check phi instead of construct distance
  // vec
  std::vector<std::vector<int>> processedCells;
  processedCells.resize(_resolution[0]);
  for (auto &row : processedCells) {
    row.resize(_resolution[1] + 1);
  }

  // Find first wavefront of surface fluid cells
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array2i ij = idToij(_id);
    int i, j;
    i = ij[0];
    j = ij[1];
    // Check if is a air surface cell
    if (_material[i][j] == Material::FluidMaterial::FLUID) {
      processedCells[i][j] = 0;
      // Look for air neighbors
      // Add air cells to processing queue
      if (i > 0 && _material[i - 1][j] == Material::FluidMaterial::AIR) {
        processedCells[i - 1][j] = 1;
        processingCells.push(ijToId(i - 1, j));
      }
      if (i < _resolution.x() - 1 &&
          _material[i + 1][j] == Material::FluidMaterial::AIR) {
        processedCells[i + 1][j] = 1;
        processingCells.push(ijToId(i + 1, j));
      }
      if (j > 0 && _material[i][j - 1] == Material::FluidMaterial::AIR) {
        processedCells[i][j - 1] = 1;
        processingCells.push(ijToId(i, j - 1));
      }
      if (j < _resolution[1] - 1 &&
          _material[i][j + 1] == Material::FluidMaterial::AIR) {
        processedCells[i][j + 1] = 1;
        processingCells.push(ijToId(i, j + 1));
      }
    }
  }

  // For each air cell, find its nearest neighbors from surface and take their
  // velocities values. In case when more than one is found, take the average
  while (!processingCells.empty()) {
    int currCell = processingCells.front();
    processingCells.pop();
    Eigen::Array2i ij = idToij(currCell);
    int i, j;
    i = ij[0];
    j = ij[1];

    if (_material[i][j] == Material::FluidMaterial::AIR &&
        _hasOppositeNeighborsWithMaterial(currCell,
                                          Material::FluidMaterial::FLUID))
      continue;

    // if (processedCells[i][j] > 20)
    // continue;

    // Find the least distance
    int leastDistace = 1e6;
    std::vector<Eigen::Array2i> neighborCells;
    if (i > 0 && processedCells[i - 1][j] <= leastDistace) {
      if (processedCells[i - 1][j] < leastDistace) {
        leastDistace = processedCells[i - 1][j];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i - 1, j);
    }
    if (i < _resolution[0] - 1 && processedCells[i + 1][j] <= leastDistace) {
      if (processedCells[i + 1][j] < leastDistace) {
        leastDistace = processedCells[i + 1][j];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i + 1, j);
    }
    if (j > 0 && processedCells[i][j - 1] <= leastDistace) {
      if (processedCells[i][j - 1] < leastDistace) {
        leastDistace = processedCells[i][j - 1];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i, j - 1);
    }
    if (j < _resolution[1] - 1 && processedCells[i][j + 1] <= leastDistace) {
      if (processedCells[i][j + 1] < leastDistace) {
        leastDistace = processedCells[i][j + 1];
        neighborCells.clear();
      }
      neighborCells.emplace_back(i, j + 1);
    }
    processedCells[i][j] = leastDistace + 1;
    // Find which faces need update
    // TODO: Change this to an enumerate
    bool faceNeedUpdate[] = {1, 1, 1, 1}; // LEFT RIGHT BOTTOM TOP
    for (auto cell : neighborCells) {
      if (cell[0] < i)
        faceNeedUpdate[0] = false;
      else if (cell[0] > i)
        faceNeedUpdate[1] = false;
      if (cell[1] < j)
        faceNeedUpdate[2] = false;
      else if (cell[1] > j)
        faceNeedUpdate[3] = false;
    }
    // Update the faces
    int nCells = neighborCells.size();
    if (faceNeedUpdate[0]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _u[cell[0]][cell[1]];
      }
      _u[i][j] = newVelocity / nCells;
      if (i > 0 && processedCells[i - 1][j] > 1e6) {
        processingCells.push(ijToId(i - 1, j));
        processedCells[i - 1][j] = processedCells[i][j] + 1;
      }
    }
    if (faceNeedUpdate[1]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _u[cell[0] + 1][cell[1]];
      }
      _u[i + 1][j] = newVelocity / nCells;
      if (i < _resolution[0] - 1 && processedCells[i + 1][j] > 1e6) {
        processingCells.push(ijToId(i + 1, j));
        processedCells[i + 1][j] = processedCells[i][j] + 1;
      }
    }
    if (faceNeedUpdate[2]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _v[cell[0]][cell[1]];
      }
      _v[i][j] = newVelocity / nCells;
      if (j > 0 && processedCells[i][j - 1] > 1e6) {
        processingCells.push(ijToId(i, j - 1));
        processedCells[i][j - 1] = processedCells[i][j] + 1;
      }
    }
    if (faceNeedUpdate[3]) {
      double newVelocity = 0.0;
      for (auto cell : neighborCells) {
        newVelocity += _v[cell[0]][cell[1] + 1];
      }
      _v[i][j + 1] = newVelocity / nCells;
      if (j < _resolution[1] - 1 && processedCells[i][j + 1] > 1e6) {
        processingCells.push(ijToId(i, j + 1));
        processedCells[i][j + 1] = processedCells[i][j] + 1;
      }
    }
  }
}

void RegularGrid2::solvePressure() {

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
    Eigen::Array2i ij = idToij(_id);
    int i, j;
    i = ij[0];
    j = ij[1];

    if (_material[i][j] != Material::FluidMaterial::FLUID)
      continue;
    int id = _getMapId(_id);

    int validCells = 0;
    std::vector<Eigen::Triplet<double>> threadTriplet;

    if (i > 0) {
      validCells++;
      if (_material[i - 1][j] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijToId(i - 1, j)),
                                   -1 / (_h.x() * _h.x()));
    }
    if (i < _resolution.x() - 1) {
      validCells++;
      if (_material[i + 1][j] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijToId(i + 1, j)),
                                   -1 / (_h.x() * _h.x()));
    }
    if (j > 0) {
      validCells++;
      if (_material[i][j - 1] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijToId(i, j - 1)),
                                   -1 / (_h.x() * _h.x()));
    }
    if (j < _resolution[1] - 1) {
      validCells++;
      if (_material[i][j + 1] != Material::FluidMaterial::AIR)
        threadTriplet.emplace_back(id, _getMapId(ijToId(i, j + 1)),
                                   -1 / (_h.x() * _h.x()));
    }

    threadTriplet.emplace_back(id, id, validCells / (_h.x() * _h.x()));
#pragma omp critical
    {
      triplets.insert(triplets.end(), threadTriplet.begin(),
                      threadTriplet.end());
    }

    divergent[id] = 0;
    divergent[id] -= (_u[i + 1][j] - _u[i][j]) / _h.x();
    divergent[id] -= (_v[i][j + 1] - _v[i][j]) / _h[1];
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
  for (int j = 0; j < _resolution[1]; j++) {
    for (int i = 1; i < _resolution.x(); i++) {
      double diffs =
          (pressure[ijToId(i, j)] - pressure[ijToId(i - 1, j)]) / _h.x();
      _u[i][j] -= _dt * diffs;
    }
  }
#pragma omp parallel for
  for (int j = 1; j < _resolution[1]; j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      double diffs =
          (pressure[ijToId(i, j)] - pressure[ijToId(i, j - 1)]) / _h[1];
      _v[i][j] -= _dt * diffs;
    }
  }
  timer.registerTime("Gradient");
  timer.evaluateComponentsTime();
}

bool RegularGrid2::_hasOppositeNeighborsWithMaterial(
    int cellId, Material::FluidMaterial material) {
  Eigen::Array2i ij = idToij(cellId);
  int i, j;
  i = ij[0];
  j = ij[1];
  if (i > 0 && _material[i - 1][j] == material && i < _resolution[0] - 2 &&
      _material[i + 1][j] == material)
    return true;
  if (j > 0 && _material[i][j - 1] == material && j < _resolution[1] - 2 &&
      _material[i][j + 1] == material)
    return true;

  return false;
}

double RegularGrid2::_interpolateVelocityU(Eigen::Array2d position) {
  double _min, _max;
  _interpolateVelocityU(position, _min, _max);
}

double RegularGrid2::_interpolateVelocityU(Eigen::Array2d position,
                                           double &_min, double &_max) {
  Eigen::Array2d h(_h[0], _h[1]);
  Eigen::Array2i resolution(_resolution[0], _resolution[1]);
  Eigen::Array2i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array2d cellCenter = index.cast<double>() * h + h / 2.0;
  // Treating values that are outside domain

  if (position[1] < cellCenter[1])
    index[1]--;

  std::vector<int> iCandidates, jCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
  }
  std::vector<Eigen::Array2d> points;
  std::vector<double> values;

  for (auto v : jCandidates)
    for (auto u : iCandidates) {
      int x, y;
      points.emplace_back(u * h[0], (v + 0.5) * h[1]);
      x = std::max(0, std::min(resolution[0], u));
      y = std::max(0, std::min(resolution[1] - 1, v));
      values.emplace_back(_u[x][y]);
      _min = std::min(values.back(), _min);
      _max = std::max(values.back(), _max);
    }
  double arrayPos[] = {position[0], position[1]};
  double velocity = Interpolator::bicubic(arrayPos, points, values);
  return velocity;
}

double RegularGrid2::_interpolateVelocityV(Eigen::Array2d position) {
  double _min, _max;
  _interpolateVelocityV(position, _min, _max);
}

double RegularGrid2::_interpolateVelocityV(Eigen::Array2d position,
                                           double &_min, double &_max) {
  Eigen::Array2d h(_h[0], _h[1]);
  Eigen::Array2i resolution(_resolution[0], _resolution[1]);
  Eigen::Array2i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array2d cellCenter = index.cast<double>() * h + h / 2.0;
  if (position[0] < cellCenter[0])
    index[0]--;

  std::vector<int> iCandidates, jCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
  }
  std::vector<Eigen::Array2d> points;
  std::vector<double> values;

  for (auto v : jCandidates)
    for (auto u : iCandidates) {
      int x, y;
      points.emplace_back((u + 0.5) * h[0], v * h[1]);
      x = std::max(0, std::min(resolution[0] - 1, u));
      y = std::max(0, std::min(resolution[1], v));
      values.emplace_back(_v[x][y]);
      _min = std::min(values.back(), _min);
      _max = std::max(values.back(), _max);
    }
  double arrayPos[] = {position[0], position[1]};
  double velocity = Interpolator::bicubic(arrayPos, points, values);
  return velocity;
}

} // namespace Ramuh
