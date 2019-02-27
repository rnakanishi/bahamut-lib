#include <structures/grid.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <omp.h>

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

void RegularGrid::advectGridVelocity() {
  Vector3d utemp[_resolution.x() + 1][_resolution.y()][_resolution.z()];
  Vector3d vtemp[_resolution.x()][_resolution.y() + 1][_resolution.z()];
  Vector3d wtemp[_resolution.x()][_resolution.y()][_resolution.z() + 1];

// Interpolate other components
#pragma omp parallel for
  // x faces - all velocities components
  for (int i = 0; i < _resolution.x() + 1; i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        double vVel = 0.0, wVel = 0.0;
        int vCount, wCount;

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
        int uCount, wCount;

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
        int uCount, vCount;

        if (k > 0) {
          uCount += 2;
          vCount += 2;
          uVel += _u[i][j][k - 1].x() + _u[i + 1][j][k - 1].x();
          vVel += _w[i][j][k - 1].y() + _w[i][j + 1][k - 1].y();
        }
        if (k < _resolution.z()) {
          uCount += 2;
          vCount += 2;
          uVel += _u[i][j][k].x() + _u[i + 1][j][k].x();
          vVel += _w[i][j][k].y() + _w[i][j + 1][k].y();
        }
        wtemp[i][j][k].x(uVel / uCount);
        wtemp[i][j][k].y(vVel / vCount);
        wtemp[i][j][k].z(_w[i][j][k].z());
      }

// Solve convection term (material derivative) - u velocities
#pragma omp parallel for
  for (int i = 1; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Vector3d backPosition,
            velPosition =
                Vector3d(i, j, k) * _h + Vector3d(0, _h.y() / 2, _h.z() / 2);
        backPosition = velPosition - utemp[i][j][k] * _dt;

        // Velocity Linear interpolation
        Vector3i index;
        double newVel = 0.0;
        int count = 0;
        index.x(std::floor(backPosition.x() / _h.x()));
        index.y(std::floor(backPosition.y() / _h.y()));
        index.z(std::floor(backPosition.z() / _h.z()));
        if (index.y() >= 0 || index.y() < _resolution.y()) {
          if (index.y() >= 0) {
            count += 2;
            newVel += utemp[i][j][k].x();
            newVel += utemp[i + 1][j][k].x();
          }
          if (index.y() < _resolution.y()) {
            count += 2;
            newVel += utemp[i][j + 1][k].x();
            newVel += utemp[i][j + 1][k].x();
          }
        }
        // TODO: Add verification for z axis
      }
// Solve convection term (material derivative) - v velocities
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 1; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Vector3d backPosition,
            velPosition =
                Vector3d(i, j, k) * _h + Vector3d(0, _h.y() / 2, _h.z() / 2);
        backPosition = velPosition - utemp[i][j][k] * _dt;

        // Velocity Linear interpolation
        Vector3i index;
        double newVel = 0.0;
        int count = 0;
        index.x(std::floor(backPosition.x() / _h.x()));
        index.y(std::floor(backPosition.y() / _h.y()));
        index.z(std::floor(backPosition.z() / _h.z()));
        if (index.x() >= 0 || index.x() < _resolution.x()) {
          if (index.x() >= 0) {
            count += 2;
            newVel += utemp[i][j][k].x();
            newVel += utemp[i][j + 1][k].x();
          }
          if (index.x() < _resolution.x()) {
            count += 2;
            newVel += utemp[i + 1][j + 1][k].x();
            newVel += utemp[i + 1][j + 1][k].x();
          }
        }
      }

  // TODO: add convective term for z axis
}

void RegularGrid::setVelocity() {

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

void RegularGrid::addGravity() {
  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = 0; j < _resolution.y() + 1; j++) {
      for (int i = 0; i < _resolution.x(); i++) {
        _v[i][j][k].y(_v[i][j][k].y() - 9.81 * _dt);
      }
    }
  }
}

void RegularGrid::solvePressure() {

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
  std::cerr << "New velocity u " << _resolution << "\n";
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
