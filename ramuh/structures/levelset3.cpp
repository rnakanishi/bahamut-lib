#include <structures/levelset3.h>
#include <geometry/matrix.h>
#include <geometry/geometry_utils.h>
#include <utils/macros.h>
#include <omp.h>
#include <queue>
#include <cmath>
#include <set>
#include <utility>
#include <algorithm>
#include <Eigen/Sparse>

namespace Ramuh {

LevelSet3::LevelSet3() : RegularGrid3() {
  _phi.resize(_resolution.x());
  for (auto &row : _phi) {
    row.resize(_resolution.y());
    for (auto &depth : row)
      depth.resize(_resolution.z(), 1e6);
  }
  _gradPhi.resize(_resolution.x() + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
  }
}

void LevelSet3::setResolution(Vector3i newResolution) {
  RegularGrid3::setResolution(newResolution);
  _phi.resize(_resolution.x());
  for (auto &row : _phi) {
    row.resize(_resolution.y());
    for (auto &depth : row)
      depth.resize(_resolution.z(), 1e6);
  }
  _gradPhi.resize(_resolution.x() + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
  }
  _velocity.resize(_resolution.x() + 1);
  for (auto &row : _velocity) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
  }
}

void LevelSet3::printVertexVelocity() {
  std::cerr << "==== x component: \n";
  for (int j = 0; j < _resolution.y() + 1; j++) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _velocity[i][j][0].x() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "==== y component: \n";
  for (int j = 0; j < _resolution.y() + 1; j++) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _velocity[i][j][0].y() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
}

void LevelSet3::addSphereSurface(Vector3d center, double radius) {
  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int i = 0; i < _resolution.x(); i++) {
        Vector3d position(Vector3i(i, j, k) * _h + _h / 2);
        double distance = (position - center).length() - radius;
        _phi[i][j][k] = distance;
      }
  checkCellMaterial();
}

void LevelSet3::addCubeSurface(Vector3d lower, Vector3d upper) {
  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int i = 0; i < _resolution.x(); i++) {
        Vector3d position(Vector3i(i, j, k) * _h + _h / 2);
        Vector3d distanceLower = (position - lower).abs();
        Vector3d distanceUpper = (position - upper).abs();
        double phi;
        if (position > lower && position < upper) {
          phi = std::min(std::fabs(_phi[i][j][k]),
                         std::min(distanceLower.min(), distanceUpper.min()));
          _phi[i][j][k] = -std::min(std::fabs(_phi[i][j][k]), phi);
        } else {
          phi = std::min(_phi[i][j][k],
                         std::min(distanceLower.min(), distanceUpper.min()));
          _phi[i][j][k] = std::min(std::fabs(_phi[i][j][k]), phi);
        }
      }
  checkCellMaterial();
}

void LevelSet3::checkCellMaterial() {
  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int i = 0; i < _resolution.x(); i++) {
        if (_phi[i][j][k] <= 0)
          _material[i][j][k] = Material::FluidMaterial::FLUID;
        else
          _material[i][j][k] = Material::FluidMaterial::AIR;
      }
}

void LevelSet3::addImplicitFunction() { NOT_IMPLEMENTED(); }

void LevelSet3::integrateLevelSet() {
  // std::vector<std::vector<double>> oldPhi;
  Matrix3<double> oldPhi;
  oldPhi.changeSize(_resolution);

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Vector3d position, backPosition, velocity, h(_h.x(), _h.y(), _h.z()),
            cellCenter;
        Vector3i index;
        double newPhi = _phi[i][j][k];
        double distanceCount = 0.0, distance = 0.0;

        position = Vector3d(i, j, k) * h + h / 2.0;
        velocity.x((_u[i][j][k].x() + _u[i + 1][j][k].x()) / 2.0);
        velocity.y((_v[i][j][k].y() + _v[i][j + 1][k].y()) / 2.0);
        velocity.z((_w[i][j][k].z() + _w[i][j][k + 1].z()) / 2.0);
        backPosition = position - (velocity * _dt);
        index.x(std::floor(backPosition.x() / h.x()));
        index.y(std::floor(backPosition.y() / h.y()));
        index.z(std::floor(backPosition.z() / h.z()));

        // Check if inside domain
        if (velocity.length() > 1e-8 && index >= Vector3i(0) &&
            index < _resolution) {
          newPhi = _interpolatePhi(Eigen::Vector3d(
              backPosition[0], backPosition[1], backPosition[2]));
        }
        oldPhi[i][j][k] = newPhi;
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++)
        _phi[i][j][k] = oldPhi[i][j][k];
}

void LevelSet3::macCormackAdvection() {
  Matrix3<double> newPhi;
  newPhi.changeSize(_resolution);
  Eigen::Array3i resolution(_resolution[0], _resolution[1], _resolution[2]);
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
        Eigen::Vector3d velocity;
        Eigen::Array3i index;
        // Find the point as if it was in previous timestep and work with it
        position = Eigen::Array3d(i, j, k).cwiseProduct(h) + (h / 2.0);
        velocity[0] = (_u[i][j][k].x() + _u[i + 1][j][k].x()) / 2.0;
        velocity[1] = (_v[i][j][k].y() + _v[i][j + 1][k].y()) / 2.0;
        velocity[2] = (_w[i][j][k].z() + _w[i][j][k + 1].z()) / 2.0;
        position -= velocity.array() * _dt;
        index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
        // TODO: look for a better solution
        if ((velocity.norm() > 1e-8) && index[0] >= 0 && index[1] >= 0 &&
            index[2] >= 0 && index[0] < _resolution[0] &&
            index[1] < _resolution[1] && index[2] < _resolution[2]) {
          double phi_n;
          phi_n = _interpolatePhi(position);
          // Advect forward in time the value there
          velocity[0] = _interpolateVelocityU(position);
          velocity[1] = _interpolateVelocityV(position);
          velocity[2] = _interpolateVelocityW(position);
          position += velocity.array() * _dt;
          double phi_n1_hat = _interpolatePhi(position);
          // Analyse the error from original to the forward advected
          double error = 0.5 * (_phi[i][j][k] - phi_n1_hat);
          newPhi[i][j][k] = phi_n + error;
        } else
          newPhi[i][j][k] = _phi[i][j][k];
      }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++)
        _phi[i][j][k] = newPhi[i][j][k];
}

double LevelSet3::_interpolatePhi(Eigen::Array3d position) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution.x(), _resolution.y(), _resolution.z());
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>().cwiseProduct(h) + h / 2.0;
  std::vector<int> iCandidates, jCandidates, kCandidates;
  iCandidates.push_back(index[0]);
  jCandidates.push_back(index[1]);
  kCandidates.push_back(index[2]);
  if (position[0] > cellCenter[0] && index[0] < resolution[0] - 1)
    iCandidates.push_back(index[0] + 1);
  else if (position[0] < cellCenter[0] && index[0] > 0)
    iCandidates.push_back(index[0] - 1);
  if (position[1] > cellCenter[1] && index[1] < resolution[1] - 1)
    jCandidates.push_back(index[1] + 1);
  else if (position[1] < cellCenter[1] && index[1] > 0)
    jCandidates.push_back(index[1] - 1);
  if (position[2] > cellCenter[2] && index[2] < resolution[2] - 1)
    kCandidates.push_back(index[2] + 1);
  else if (position[2] < cellCenter[2] && index[2] > 0)
    kCandidates.push_back(index[2] - 1);

  // Catmull-Rom like interpolation
  double newPhi = 0., distance = 0., distanceCount = 0.;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        Eigen::Array3d centerPosition = Eigen::Array3d(u, v, w) * h + h / 2.0;
        if (position.matrix().isApprox(centerPosition.matrix(), 1e-8))
          return _phi[u][v][w];
        distance = (position - centerPosition).matrix().norm();
        distanceCount += 1. / distance;
        newPhi += _phi[u][v][w] / distance;
      }
  return newPhi / distanceCount;
}

std::vector<std::vector<double>> &LevelSet3::operator[](const int i) {
  return _phi[i];
}

void LevelSet3::solvePressure() {
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
        double h2 = _h.x() * _h.x();
        double centerWeight = 0.0;
        int id = ijkToId(i, j, k);
        std::vector<Eigen::Triplet<double>> threadTriplet;
        Eigen::Array3d center(i * _h[0], j * _h[1], k * _h[2]);

        if (i > 0) {
          validCells++;
          if (_material[i - 1][j][k] != Material::FluidMaterial::AIR) {
            centerWeight += 1 / h2;
            threadTriplet.emplace_back(id, ijkToId(i - 1, j, k), -1 / (h2));
          } else {
            double theta;
            Eigen::Array3d surface = _findSurfaceCoordinate(
                Eigen::Array3i(i, j, k), Eigen::Array3i(i - 1, j, k));
            theta = (surface - center).matrix().norm() / _h.x();
            centerWeight += 1 / (h2 * theta);
          }
        }
        if (i < _resolution.x() - 1) {
          validCells++;
          if (_material[i + 1][j][k] != Material::FluidMaterial::AIR) {
            centerWeight += 1 / h2;
            threadTriplet.emplace_back(id, ijkToId(i + 1, j, k), -1 / (h2));
          } else {
            double theta;
            Eigen::Array3d surface = _findSurfaceCoordinate(
                Eigen::Array3i(i, j, k), Eigen::Array3i(i + 1, j, k));
            theta = (surface - center).matrix().norm() / _h.x();
            centerWeight += 1 / (h2 * theta);
          }
        }
        if (j > 0) {
          validCells++;
          if (_material[i][j - 1][k] != Material::FluidMaterial::AIR) {
            centerWeight += 1 / h2;
            threadTriplet.emplace_back(id, ijkToId(i, j - 1, k), -1 / (h2));
          } else {
            double theta;
            Eigen::Array3d surface = _findSurfaceCoordinate(
                Eigen::Array3i(i, j, k), Eigen::Array3i(i, j - 1, k));
            theta = (surface - center).matrix().norm() / _h.x();
            centerWeight += 1 / (h2 * theta);
          }
        }
        if (j < _resolution.y() - 1) {
          validCells++;
          if (_material[i][j + 1][k] != Material::FluidMaterial::AIR) {
            centerWeight += 1 / h2;
            threadTriplet.emplace_back(id, ijkToId(i, j + 1, k), -1 / (h2));
          } else {
            double theta;
            Eigen::Array3d surface = _findSurfaceCoordinate(
                Eigen::Array3i(i, j, k), Eigen::Array3i(i, j + 1, k));
            theta = (surface - center).matrix().norm() / _h.x();
            centerWeight += 1 / (h2 * theta);
          }
        }
        if (k > 0) {
          validCells++;
          if (_material[i][j][k - 1] != Material::FluidMaterial::AIR) {
            centerWeight += 1 / h2;
            threadTriplet.emplace_back(id, ijkToId(i, j, k - 1), -1 / (h2));
          } else {
            double theta;
            Eigen::Array3d surface = _findSurfaceCoordinate(
                Eigen::Array3i(i, j, k), Eigen::Array3i(i, j, k - 1));
            theta = (surface - center).matrix().norm() / _h.x();
            centerWeight += 1 / (h2 * theta);
          }
        }
        if (k < _resolution.z() - 1) {
          validCells++;
          if (_material[i][j][k + 1] != Material::FluidMaterial::AIR) {
            centerWeight += 1 / h2;
            threadTriplet.emplace_back(id, ijkToId(i, j, k + 1), -1 / (h2));
          } else {
            double theta;
            Eigen::Array3d surface = _findSurfaceCoordinate(
                Eigen::Array3i(i, j, k), Eigen::Array3i(i, j, k + 1));
            theta = (surface - center).matrix().norm() / _h.x();
            centerWeight += 1 / (h2 * theta);
          }
        }

        threadTriplet.emplace_back(id, id, validCells / (h2));
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
        _u[i][j][k].x(_u[i][j][k].x() - _dt *
                                            (pressure[ijkToId(i, j, k)] -
                                             pressure[ijkToId(i - 1, j, k)]) /
                                            _h.x());
      }
    }
  }
#pragma omp parallel for
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
#pragma omp parallel for
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

double LevelSet3::_solveEikonal(glm::ivec3 cellId) {
  int i, j, k;
  i = cellId[0];
  j = cellId[1];
  k = cellId[2];
  std::vector<double> distances(3, 0);
  distances[0] =
      std::min(std::fabs(_phi[std::max(0, i - 1)][j][k]),
               std::fabs(_phi[std::min(_resolution.x() - 1, i + 1)][j][k]));
  distances[1] =
      std::min(std::fabs(_phi[i][std::max(0, j - 1)][k]),
               std::fabs(_phi[i][std::min(_resolution.y() - 1, j + 1)][k]));
  distances[2] =
      std::min(std::fabs(_phi[i][j][std::max(0, k - 1)]),
               std::fabs(_phi[i][j][std::min(_resolution.z() - 1, k + 1)]));

  std::sort(distances.begin(), distances.end());

  double newPhi = distances[0] + _h.x();
  if (newPhi > distances[1]) {
    double h2 = _h.x() * _h.x();
    newPhi = 0.5 * (distances[0] + distances[1] +
                    std::sqrt(2 * h2 - ((distances[1] - distances[0]) *
                                        (distances[1] - distances[0]))));
    if (std::fabs(newPhi) > distances[2]) {
      newPhi = (distances[0] + distances[1] + distances[2]);
      newPhi += std::sqrt(
          std::max(0.0, (distances[0] + distances[1] + distances[2]) *
                                (distances[0] + distances[1] + distances[2]) -
                            3 * (distances[0] * distances[0] +
                                 distances[1] * distances[1] +
                                 distances[2] * distances[2] - h2)));
      newPhi /= 3;
    }
  }
  if (std::fabs(newPhi) < 1e-8 || std::isnan(newPhi) || newPhi == 0 ||
      newPhi > 5)
    float np = newPhi;
  return newPhi;
}

void LevelSet3::redistance() {
  Matrix3<double> tempPhi;
  Matrix3<bool> processed;
  std::priority_queue<std::pair<double, int>,
                      std::vector<std::pair<double, int>>,
                      std::greater<std::pair<double, int>>>
      cellsQueue;
  std::set<int> cellsAdded;

  tempPhi.changeSize(_resolution, 1e8);
  processed.changeSize(_resolution, false);

#pragma omp parallel for
  // Find surface cells adding them to queue
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        glm::vec3 position = glm::vec3(i * _h.x(), j * _h.y(), k * _h.z());
        glm::vec3 intersections[3];
        int nintersecs = 0;
        double cellPhi = _phi[i][j][k];
        int cellSign = cellPhi / std::fabs(cellPhi);
        if (cellPhi == 0)
          cellSign = 1;

        bool isSurface = false, intersected = false;
        double theta = 1e8;
        if (i < _resolution.x() - 1 &&
            std::signbit(cellSign) != std::signbit(_phi[i + 1][j][k])) {
          isSurface = true;
          intersected = true;
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i + 1][j][k])));
          intersections[nintersecs] =
              glm::vec3(position[0] + theta * _h.x(), position[1], position[2]);
        }
        if (i > 0 &&
            std::signbit(cellSign) != std::signbit(_phi[i - 1][j][k])) {
          isSurface = true;
          intersected = true;
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i - 1][j][k])));
          intersections[nintersecs] =
              glm::vec3(position[0] - theta * _h.x(), position[1], position[2]);
        }
        theta = 1e8;
        if (intersected) {
          intersected = false;
          nintersecs++;
        }
        if (j < _resolution.y() - 1 &&
            std::signbit(cellSign) != std::signbit(_phi[i][j + 1][k])) {
          isSurface = true;
          intersected = true;
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i][j + 1][k])));
          intersections[nintersecs] =
              glm::vec3(position[0], position[1] + theta * _h.y(), position[2]);
        }
        if (j > 0 &&
            std::signbit(cellSign) != std::signbit(_phi[i][j - 1][k])) {
          isSurface = true;
          intersected = true;
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i][j - 1][k])));
          intersections[nintersecs] =
              glm::vec3(position[0], position[1] - theta * _h.y(), position[2]);
        }
        theta = 1e8;
        if (intersected) {
          intersected = false;
          nintersecs++;
        }
        if (k < _resolution.z() - 1 &&
            std::signbit(cellSign) != std::signbit(_phi[i][j][k + 1])) {
          isSurface = true;
          intersected = true;
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i][j][k + 1])));
          intersections[nintersecs] =
              glm::vec3(position[0], position[1], position[2] + theta * _h.z());
        }
        if (k > 0 &&
            std::signbit(cellSign) != std::signbit(_phi[i][j][k - 1])) {
          isSurface = true;
          intersected = true;
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i][j][k - 1])));
          intersections[nintersecs] =
              glm::vec3(position[0], position[1], position[2] - theta * _h.z());
        }
        if (intersected) {
          intersected = false;
          nintersecs++;
        }

        if (isSurface) {
          glm::vec3 proj;
          if (nintersecs == 1) {
            proj = intersections[0];
          } else if (nintersecs == 2) {
            proj = Geometry::closestPointPlane(position, intersections[0],
                                               intersections[1]);
          } else if (nintersecs == 3) {
            proj = Geometry::closestPointTriangle(
                position, intersections[0], intersections[1], intersections[2]);
          }
          double distance = glm::length(proj - position);
          if (std::isnan(distance) || distance == 0 || distance > 5)
            float d = distance;
          if (std::fabs(distance) < 1e-8)
            distance = 0;

          tempPhi[i][j][k] = cellSign * distance;
          // _phi[i][j][k];
          // cellSign * _solveEikonal(glm::ivec3(i, j, k));
          // cellSign *std::min(std::fabs(tempPhi[i][j][k]), theta * _h.x());
          processed[i][j][k] = true;
#pragma omp critical
          {
            cellsQueue.push(
                std::make_pair(std::fabs(tempPhi[i][j][k]), ijkToId(i, j, k)));
            cellsAdded.insert(ijkToId(i, j, k));
          }
        }
      }

  // Propagating distances
  while (!cellsQueue.empty()) {
    int cellId = cellsQueue.top().second;
    cellsQueue.pop();
    int k = cellId / (_resolution.x() * _resolution.y());
    int j = (cellId % (_resolution.x() * _resolution.y())) / _resolution.x();
    int i = (cellId % (_resolution.x() * _resolution.y())) % _resolution.x();
    // TODO: check this caculation

    if (!processed[i][j][k]) {
      processed[i][j][k] = true;
      int distanceSignal = _phi[i][j][k] / std::fabs(_phi[i][j][k]);
      double newPhi = _solveEikonal(glm::ivec3(i, j, k));

      if (newPhi < std::fabs(tempPhi[i][j][k]))
        tempPhi[i][j][k] = distanceSignal * newPhi;
    }
    // Add all neighbors of the cell to queue
    if (i > 0 && cellsAdded.find(ijkToId(i - 1, j, k)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i - 1][j][k]), ijkToId(i - 1, j, k)));
      cellsAdded.insert(ijkToId(i - 1, j, k));
    }
    if (i < _resolution.x() - 1 &&
        cellsAdded.find(ijkToId(i + 1, j, k)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i + 1][j][k]), ijkToId(i + 1, j, k)));
      cellsAdded.insert(ijkToId(i + 1, j, k));
    }
    if (j > 0 && cellsAdded.find(ijkToId(i, j - 1, k)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i][j - 1][k]), ijkToId(i, j - 1, k)));
      cellsAdded.insert(ijkToId(i, j - 1, k));
    }
    if (j < _resolution.y() - 1 &&
        cellsAdded.find(ijkToId(i, j + 1, k)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i][j + 1][k]), ijkToId(i, j + 1, k)));
      cellsAdded.insert(ijkToId(i, j + 1, k));
    }
    if (k > 0 && cellsAdded.find(ijkToId(i, j, k - 1)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i][j][k - 1]), ijkToId(i, j, k - 1)));
      cellsAdded.insert(ijkToId(i, j, k - 1));
    }
    if (k < _resolution.z() - 1 &&
        cellsAdded.find(ijkToId(i, j, k + 1)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i][j][k + 1]), ijkToId(i, j, k + 1)));
      cellsAdded.insert(ijkToId(i, j, k + 1));
    }
  }

  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++)
        _phi[i][j][k] = tempPhi[i][j][k];
}

void LevelSet3::printLevelSetValue() {

  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = _resolution.y() - 1; j >= 0; j--) {
      for (int i = 0; i < _resolution.x(); i++) {
        if (k == 15)
          std::cout << _phi[i][j][k] << ' ';
      }
      if (k == 15)
        std::cout << std::endl;
    }
    if (k == 15)
      std::cout << std::endl;
  }
}

TriangleMesh LevelSet3::marchingTetrahedra() {
  TriangleMesh mesh;

  // Look for octants that have different levelset signal
  for (int i = 0; i < _resolution.x() - 1; i++)
    for (int j = 0; j < _resolution.y() - 1; j++)
      for (int k = 0; k < _resolution.z() - 1; k++) {
        std::vector<glm::ivec3> tetraIndex;
        // Tetrahedron 0127
        tetraIndex.clear();
        tetraIndex.emplace_back(i, j, k);
        tetraIndex.emplace_back(i + 1, j, k);
        tetraIndex.emplace_back(i, j + 1, k);
        tetraIndex.emplace_back(i + 1, j + 1, k + 1);
        _triangulate(tetraIndex, mesh);

        // Tetrahedron 0456
        tetraIndex.clear();
        tetraIndex.emplace_back(i, j, k);
        tetraIndex.emplace_back(i, j, k + 1);
        tetraIndex.emplace_back(i + 1, j, k + 1);
        tetraIndex.emplace_back(i, j + 1, k + 1);
        _triangulate(tetraIndex, mesh);

        // Tetrahedron 0267
        tetraIndex.clear();
        tetraIndex.emplace_back(i, j, k);
        tetraIndex.emplace_back(i, j + 1, k);
        tetraIndex.emplace_back(i, j + 1, k + 1);
        tetraIndex.emplace_back(i + 1, j + 1, k + 1);
        _triangulate(tetraIndex, mesh);

        // Tetrahedron 0157
        tetraIndex.clear();
        tetraIndex.emplace_back(i, j, k);
        tetraIndex.emplace_back(i + 1, j, k);
        tetraIndex.emplace_back(i + 1, j, k + 1);
        tetraIndex.emplace_back(i + 1, j + 1, k + 1);
        _triangulate(tetraIndex, mesh);

        // Tetrahedron 1237
        tetraIndex.clear();
        tetraIndex.emplace_back(i + 1, j, k);
        tetraIndex.emplace_back(i, j + 1, k);
        tetraIndex.emplace_back(i + 1, j + 1, k);
        tetraIndex.emplace_back(i + 1, j + 1, k + 1);
        _triangulate(tetraIndex, mesh);

        // Tetrahedron 0567
        tetraIndex.clear();
        tetraIndex.emplace_back(i, j, k);
        tetraIndex.emplace_back(i + 1, j, k + 1);
        tetraIndex.emplace_back(i, j + 1, k + 1);
        tetraIndex.emplace_back(i + 1, j + 1, k + 1);
        _triangulate(tetraIndex, mesh);
      }
  // Prepare tetrahedra decomposition and triangulate them

  return mesh;
}

void LevelSet3::_triangulate(std::vector<glm::ivec3> vertices,
                             TriangleMesh &mesh) {
  int triIndex = 0;
  if (_phi[vertices[0][0]][vertices[0][1]][vertices[0][2]] < 0)
    triIndex |= 1;
  if (_phi[vertices[1][0]][vertices[1][1]][vertices[1][2]] < 0)
    triIndex |= 2;
  if (_phi[vertices[2][0]][vertices[2][1]][vertices[2][2]] < 0)
    triIndex |= 4;
  if (_phi[vertices[3][0]][vertices[3][1]][vertices[3][2]] < 0)
    triIndex |= 8;

  glm::vec3 vertex;
  glm::ivec3 face;
  // 15 total cases to check
  // TODO: Check normal direction and change vertices order accordingly
  switch (triIndex) {
  case 0:
  case 15:
    break;
  case 1:
  case 14:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[1]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[2]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[3]));
    mesh.addFace(face);
    break;
  case 2:
  case 13:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[2]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[3]));
    mesh.addFace(face);
    break;
  case 4:
  case 11:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[1]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[3]));
    mesh.addFace(face);
    break;
  case 8:
  case 7:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[1]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[2]));
    mesh.addFace(face);
    break;
  case 3:
  case 12:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[2]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[3]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[3]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[2]));
    mesh.addFace(face);
    break;
  case 5:
  case 10:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[1]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[3]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[3]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[2]));
    mesh.addFace(face);
    break;
  case 9:
  case 6:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[1]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[2]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[3]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[3]));
    mesh.addFace(face);
    break;
  }
}

// TODO: remove this function
glm::vec3 LevelSet3::_findSurfaceCoordinate(glm::ivec3 v1, glm::ivec3 v2) {
  glm::vec3 coord1, coord2, direction;
  coord1 = glm::vec3(v1[0] * _h.x(), v1[1] * _h.y(), v1[2] * _h.z());
  coord2 = glm::vec3(v2[0] * _h.x(), v2[1] * _h.y(), v2[2] * _h.z());
  direction = coord2 - coord1;

  float phi1, phi2;
  phi1 = _phi[v1[0]][v1[1]][v1[2]];
  phi2 = _phi[v2[0]][v2[1]][v2[2]];
  return coord1 - phi1 * direction / (phi2 - phi1);
}

Eigen::Array3d LevelSet3::_findSurfaceCoordinate(Eigen::Array3i v1,
                                                 Eigen::Array3i v2) {
  Eigen::Array3d coord1, coord2, direction;
  coord1 = Eigen::Array3d(v1[0] * _h.x(), v1[1] * _h.y(), v1[2] * _h.z());
  coord2 = Eigen::Array3d(v2[0] * _h.x(), v2[1] * _h.y(), v2[2] * _h.z());
  direction = coord2 - coord1;

  float phi1, phi2;
  phi1 = _phi[v1[0]][v1[1]][v1[2]];
  phi2 = _phi[v2[0]][v2[1]][v2[2]];
  return coord1 - phi1 * direction / (phi2 - phi1);
}

} // namespace Ramuh