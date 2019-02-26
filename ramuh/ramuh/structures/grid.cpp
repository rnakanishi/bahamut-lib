#include <structures/grid.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
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

  // _pressure.resize(_resolution.x());
  // for (auto &row : _pressure) {
  //   row.resize(_resolution.y());
  //   for (auto &depth : row)
  //     depth.resize(_resolution.z());
  // }

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

void RegularGrid::solvePressure() {

  // Compute velocity divergent over cell center
  // std::vector<double> divergent(cellCount(), 0.0);
  int nCells = cellCount();
  Eigen::VectorXd divergent;
  Eigen::SparseMatrix<double> pressureMatrix(nCells, nCells);
  std::vector<Eigen::Triplet<double>> triplets;

  divergent = Eigen::VectorXd::Zero(nCells);

#pragma omp parallel
  // Solve pressure Poisson equation
  {
#pragma omp for
    for (int k = 0; k < _resolution.z(); k++) {
      for (int j = 0; j < _resolution.y(); j++) {
        for (int i = 0; i < _resolution.x(); i++) {
          int validCells = 0;
          int id = ijkToId(i, j, k);
          std::vector<Eigen::Triplet<double>> threadTriplet;
          // TODO: create a validity check for fluid cells
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

  // Correct velocity through pressure gradient
}

} // namespace Ramuh
