#include <structures/grid.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

namespace Ramuh {

RegularGrid::RegularGrid() {
  setResolution(Vector3i(32, 32, 32));
  _h = Vector3d(1. / 32);
  _domainSize = Vector3d(1.0);
}

Vector3d RegularGrid::domainSize() { return _domainSize; }

Vector3i RegularGrid::resolution() { return _resolution; }

Vector3d RegularGrid::h() { return _h; }

Vector3d RegularGrid::gridSize() { return _domainSize; }

void RegularGrid::setSize(Vector3d newSize) {
  _domainSize = newSize;
  setH(_domainSize / _resolution);
}

void RegularGrid::setResolution(Vector3i newResolution) {
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

void RegularGrid::setH(Vector3d newH) { _h = newH; }

void RegularGrid::setVelocity() {

  for (int i = 0; i < _resolution.x() + 1; i++) {
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _u[i][j][k] = 0;
      }
  }

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int k = 0; k < _resolution.z(); k++) {
        _v[i][j][k] = 0;
      }
  }

  for (int i = 0; i < _resolution.x(); i++) {
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z() + 1; k++) {
        _w[i][j][k] = 0;
      }
  }
}

int RegularGrid::cellCount() {
  return _resolution.x() * _resolution.y() * _resolution.z();
}

int RegularGrid::ijkToId(int i, int j, int k) {
  return k * _resolution.x() * _resolution.y() + j * _resolution.x() + i;
}

Vector3i RegularGrid::idToijk(int id) {}

void RegularGrid::printFaceVelocity() {

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

void RegularGrid::boundaryVelocities() {
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

void RegularGrid::addGravity() {
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y() + 1; j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        _v[i][j][k] = _v[i][j][k] - 9.81 * _dt;
      }
    }
  }
}

void RegularGrid::solvePressure() {

  // Compute velocity divergent over cell center
  // std::vector<double> divergent(cellCount(), 0.0);
  int nCells = cellCount();
  Eigen::VectorXd divergent, pressure;
  Eigen::SparseMatrix<double> pressureMatrix(nCells, nCells);
  std::vector<Eigen::Triplet<double>> triplets;

  divergent = Eigen::VectorXd::Zero(nCells);
  divergent = Eigen::VectorXd::Zero(nCells);

#pragma omp parallel
  // Solve pressure Poisson equation
  {
#pragma omp for
    for (int k = 0; k < _resolution.z(); k++) {
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
              threadTriplet.emplace_back(id, ijkToId(i - 1, j, k), 1);
            }
            if (i < _resolution.x() - 1) {
              validCells++;
              threadTriplet.emplace_back(id, ijkToId(i + 1, j, k), 1);
            }
          }
          if (j > 0 || j < _resolution.y() - 1) {
            if (j > 0) {
              validCells++;
              threadTriplet.emplace_back(id, ijkToId(i, j - 1, k), 1);
            }
            if (j < _resolution.y() - 1) {
              validCells++;
              threadTriplet.emplace_back(id, ijkToId(i, j + 1, k), 1);
            }
          }
          if (k > 0 || k < _resolution.z() - 1) {
            if (k > 0) {
              validCells++;
              threadTriplet.emplace_back(id, ijkToId(i, j, k - 1), 1);
            }
            if (k < _resolution.z() - 1) {
              validCells++;
              threadTriplet.emplace_back(id, ijkToId(i, j, k + 1), 1);
            }
          }
          threadTriplet.emplace_back(id, ijkToId(i, j, k), -validCells);
#pragma omp critical
          {
            triplets.insert(triplets.end(), threadTriplet.begin(),
                            threadTriplet.end());
          }

          divergent[ijkToId(i, j, k)] = 0;
          divergent[ijkToId(i, j, k)] +=
              (_u[i + 1][j][k] - _u[i][j][k]) / _h.x();
          divergent[ijkToId(i, j, k)] +=
              (_u[i][j + 1][k] - _u[i][j][k]) / _h.y();
          divergent[ijkToId(i, j, k)] +=
              (_u[i][j][k + 1] - _u[i][j][k]) / _h.z();
        }
      }
    }
  }

  pressureMatrix.setFromTriplets(triplets.begin(), triplets.end());

  // SOlve pressure Poisson system
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
  solver.compute(pressureMatrix);
  pressure = solver.solve(divergent);

  // Correct velocity through pressure gradient

  for (int k = 1; k < _resolution.z(); k++) {
    for (int j = 1; j < _resolution.y(); j++) {
      for (int i = 1; i < _resolution.x(); i++) {
        _u[i][j][k] = _u[i][j][k] - (pressure[ijkToId(i, j, k)] -
                                     pressure[ijkToId(i - 1, j, k)]) /
                                        _h.x();
        _v[i][j][k] = _v[i][j][k] - (pressure[ijkToId(i, j, k)] -
                                     pressure[ijkToId(i, j - 1, k)]) /
                                        _h.y();
        _w[i][j][k] = _w[i][j][k] - (pressure[ijkToId(i, j, k)] -
                                     pressure[ijkToId(i, j, k - 1)]) /
                                        _h.z();
      }
    }
  }
}

} // namespace Ramuh
