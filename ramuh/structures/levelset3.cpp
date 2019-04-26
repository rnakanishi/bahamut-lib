#include <structures/levelset3.h>
#include <geometry/matrix.h>
#include <geometry/geometry_utils.h>
#include <utils/macros.h>
#include <utils/timer.hpp>
#include <omp.h>
#include <queue>
#include <cmath>
#include <set>
#include <utility>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <iomanip>
#include <fstream>

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

  _isPressureSecondOrder = true;
}

void LevelSet3::setPressureSecondOrder(bool value) {
  _isPressureSecondOrder = value;
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
        // Find threshold values for clamping
        double clamp[2]; // 0: min, 1: max
        clamp[0] = 1e8;
        clamp[1] = -1e8;
        if (i > 0) {
          clamp[0] = std::min(clamp[0], _phi[i - 1][j][k]);
          clamp[1] = std::max(clamp[1], _phi[i - 1][j][k]);
        }
        if (i < resolution[0] - 1) {
          clamp[0] = std::min(clamp[0], _phi[i + 1][j][k]);
          clamp[1] = std::max(clamp[1], _phi[i + 1][j][k]);
        }
        if (j > 0) {
          clamp[0] = std::min(clamp[0], _phi[i][j - 1][k]);
          clamp[1] = std::max(clamp[1], _phi[i][j - 1][k]);
        }
        if (j < resolution[1] - 1) {
          clamp[0] = std::min(clamp[0], _phi[i][j + 1][k]);
          clamp[1] = std::max(clamp[1], _phi[i][j + 1][k]);
        }
        if (k > 0) {
          clamp[0] = std::min(clamp[0], _phi[i][j][k - 1]);
          clamp[1] = std::max(clamp[1], _phi[i][j][k - 1]);
        }
        if (k < resolution[2] - 1) {
          clamp[0] = std::min(clamp[0], _phi[i][j][k + 1]);
          clamp[1] = std::max(clamp[1], _phi[i][j][k + 1]);
        }
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
          newPhi[i][j][k] =
              std::max(clamp[0], std::min(clamp[1], phi_n + error));
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
  if (index[0] >= 0 && index[0] < _resolution[0])
    iCandidates.push_back(index[0]);
  if (index[1] >= 0 && index[1] < _resolution[1])
    jCandidates.push_back(index[1]);
  if (index[2] >= 0 && index[2] < _resolution[2])
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
  Eigen::Array2d clamp; // 0: min, 1: max
  clamp[0] = 1e8;
  clamp[1] = -1e8;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        Eigen::Array3d centerPosition = Eigen::Array3d(u, v, w) * h + h / 2.0;
        distance = (position - centerPosition).matrix().norm();
        if (distance < 1e-6)
          return _phi[u][v][w];
        distanceCount += 1. / distance;
        newPhi += _phi[u][v][w] / distance;
        clamp[0] = std::min(clamp[0], _phi[u][v][w]);
        clamp[1] = std::max(clamp[1], _phi[u][v][w]);
      }
  if (distanceCount == 0 || distanceCount > 1e8)
    throw("LevelSet3::interpolatePhi: distanceCount zero or infinity\n");
  return std::max(clamp[0], std::min(clamp[1], newPhi / distanceCount));
}

std::vector<std::vector<double>> &LevelSet3::operator[](const int i) {
  return _phi[i];
}

void LevelSet3::solvePressure() {
  // Compute velocity divergent over cell center
  // std::vector<double> divergent(cellCount(), 0.0);
  // TODO: Change to fluid cells count
  Eigen::VectorXd divergent, pressure;
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  std::vector<Eigen::Triplet<double>> triplets;

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

  Timer timer;
  // Solve pressure Poisson equation
  divergent = Eigen::VectorXd::Zero(cellCount());
  int ghostId = _getMapId(-1);
#pragma omp for
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
    double h2 = _h.x() * _h.x();
    double centerWeight = 0.0;
    // int id = ijkToId(i, j, k);
    std::vector<Eigen::Triplet<double>> threadTriplet;
    Eigen::Array3d center = Eigen::Array3d(i, j, k).cwiseProduct(h) + h / 2.0;
    // (i * _h[0], j * _h[1], k * _h[2]);

    if (i > 0) {
      validCells++;
      if (_material[i - 1][j][k] != Material::FluidMaterial::AIR) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i - 1, j, k)),
                                   -1 / (h2));
      } else {
        double theta;
        Eigen::Array3d surface = _findSurfaceCoordinate(
            Eigen::Array3i(i, j, k), Eigen::Array3i(i - 1, j, k));
        theta = (surface - center).matrix().norm() / _h.x();
        if (_isPressureSecondOrder) {
          centerWeight += (1 / (h2 * theta));
          threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (i < _resolution.x() - 1) {
      validCells++;
      if (_material[i + 1][j][k] != Material::FluidMaterial::AIR) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i + 1, j, k)),
                                   -1 / (h2));
      } else {
        double theta;
        Eigen::Array3d surface = _findSurfaceCoordinate(
            Eigen::Array3i(i, j, k), Eigen::Array3i(i + 1, j, k));
        theta = (surface - center).matrix().norm() / _h.x();
        if (_isPressureSecondOrder) {
          centerWeight += (1 / (h2 * theta));
          threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (j > 0) {
      validCells++;
      if (_material[i][j - 1][k] != Material::FluidMaterial::AIR) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j - 1, k)),
                                   -1 / (h2));
      } else {
        double theta;
        Eigen::Array3d surface = _findSurfaceCoordinate(
            Eigen::Array3i(i, j, k), Eigen::Array3i(i, j - 1, k));
        theta = (surface - center).matrix().norm() / _h.y();
        if (_isPressureSecondOrder) {
          centerWeight += (1 / (h2 * theta));
          threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (j < _resolution.y() - 1) {
      validCells++;
      if (_material[i][j + 1][k] != Material::FluidMaterial::AIR) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j + 1, k)),
                                   -1 / (h2));
      } else {
        double theta;
        Eigen::Array3d surface = _findSurfaceCoordinate(
            Eigen::Array3i(i, j, k), Eigen::Array3i(i, j + 1, k));
        theta = (surface - center).matrix().norm() / _h.y();
        if (_isPressureSecondOrder) {
          centerWeight += (1 / (h2 * theta));
          threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (k > 0) {
      validCells++;
      if (_material[i][j][k - 1] != Material::FluidMaterial::AIR) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j, k - 1)),
                                   -1 / (h2));
      } else {
        double theta;
        Eigen::Array3d surface = _findSurfaceCoordinate(
            Eigen::Array3i(i, j, k), Eigen::Array3i(i, j, k - 1));
        theta = (surface - center).matrix().norm() / _h.z();
        if (_isPressureSecondOrder) {
          centerWeight += (1 / (h2 * theta));
          threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (k < _resolution.z() - 1) {
      validCells++;
      if (_material[i][j][k + 1] != Material::FluidMaterial::AIR) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j, k + 1)),
                                   -1 / (h2));
      } else {
        double theta;
        Eigen::Array3d surface = _findSurfaceCoordinate(
            Eigen::Array3i(i, j, k), Eigen::Array3i(i, j, k + 1));
        theta = (surface - center).matrix().norm() / _h.z();
        if (_isPressureSecondOrder) {
          centerWeight += (1 / (h2 * theta));
          threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (std::isinf(centerWeight))
      std::cerr << ">>>>>>>>>>>>Inifinte!\n";
    if (std::isnan(centerWeight))
      std::cerr << "NaN\n";
    threadTriplet.emplace_back(id, id, centerWeight);
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
  divergent[ghostId] = 0;
  triplets.emplace_back(ghostId, ghostId, 1);

  int nCells = idMap.size();
  pressure = Eigen::VectorXd::Zero(nCells);
  // divergent.conservativeResize(nCells);
  divergent.conservativeResize(nCells);

  Eigen::SparseMatrix<double> pressureMatrix(nCells, nCells);
  // pressureMatrix.setIdentity();
  pressureMatrix.setFromTriplets(triplets.begin(), triplets.end());
  timer.registerTime("Assembly");
  if (pressureMatrix.toDense().hasNaN()) {
    std::cerr << "NaN values on matrix\n";
    exit(-12);
  }

  // SOlve pressure Poisson system
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                  Eigen::DiagonalPreconditioner<double>>
      solver;
  // solver.setMaxIterations(100);
  solver.compute(pressureMatrix);
  timer.registerTime("Prepare Matrix");
  pressure = solver.solve(divergent);
  timer.registerTime("System");

  std::cerr << "Solved with " << solver.iterations() << " iterations.\n";
  std::cerr << "Solver error: " << solver.error() << std::endl;
  {
    Eigen::VectorXd x = pressureMatrix * pressure;
    std::cerr << "L2 norm: " << (x - divergent).norm() << std::endl;
    if (!x.isApprox(divergent, 1e-2))
      std::cerr << "Error is too high\n";
  }

  std::ofstream file("eigenMatrix");
  if (file.is_open()) {
    file << pressureMatrix;
  }
  file.close();
  file.open("divergent");
  if (file.is_open()) {
    file << divergent;
  }
  file.close();
  file.open("pressure");
  if (file.is_open()) {
    file << pressure;
  }
  file.close();

  // Correct velocity through pressure gradient
  _maxVelocity[0] = _maxVelocity[1] = _maxVelocity[2] = -1e8;
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    Eigen::Array3i ijk = idToijk(id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    if (_material[i][j][k] != Material::FluidMaterial::FLUID) {
      bool hasFluidNeighbor = false;
      if (i > 0 && _material[i - 1][j][k] == Material::FluidMaterial::FLUID)
        hasFluidNeighbor = true;
      if (j > 0 && _material[i][j - 1][k] == Material::FluidMaterial::FLUID)
        hasFluidNeighbor = true;
      if (k > 0 && _material[i][j][k - 1] == Material::FluidMaterial::FLUID)
        hasFluidNeighbor = true;
      if (!hasFluidNeighbor)
        continue;
    }
    double pressure1, pressure2;
    auto it = idMap.find(id);

    pressure1 = (it == idMap.end()) ? 0.0 : pressure[it->second];
    pressure2 = 0.0;
    if (i > 0) {
      auto it = idMap.find(ijkToId(i - 1, j, k));
      if (it == idMap.end())
        pressure2 = 0.0;
      else
        pressure2 = pressure[it->second];
      _u[i][j][k].x(_u[i][j][k].x() - _dt * (pressure1 - pressure2) / _h.x());
      if (std::isnan(_u[i][j][k].x()) || std::isinf(_u[i][j][k].x())) {
        std::cerr << "Infinite velocity component";
      }
      _maxVelocity[0] = std::max(_maxVelocity[0], std::fabs(_u[i][j][k][0]));
    }
    pressure2 = 0.0;
    if (j > 0) {
      auto it = idMap.find(ijkToId(i, j - 1, k));
      if (it == idMap.end())
        pressure2 = 0.0;
      else
        pressure2 = pressure[it->second];
      _v[i][j][k].y(_v[i][j][k].y() - _dt * (pressure1 - pressure2) / _h.y());
      if (std::isnan(_v[i][j][k].y()) || std::isinf(_v[i][j][k].y())) {
        std::cerr << "Infinite velocity component";
      }
      _maxVelocity[1] = std::max(_maxVelocity[1], std::fabs(_v[i][j][k][1]));
    }
    pressure2 = 0.0;
    if (k > 0) {
      auto it = idMap.find(ijkToId(i, j, k - 1));
      if (it == idMap.end())
        pressure2 = 0.0;
      else
        pressure2 = pressure[it->second];
      _w[i][j][k].z(_w[i][j][k].z() - _dt * (pressure1 - pressure2) / _h.z());
      if (std::isnan(_w[i][j][k].z()) || std::isinf(_w[i][j][k].z())) {
        std::cerr << "Infinite velocity component";
      }
      _maxVelocity[2] = std::max(_maxVelocity[2], std::fabs(_w[i][j][k][2]));
    }
  }
  if (_maxVelocity[0] > 1e4 || _maxVelocity[1] > 1e4 || _maxVelocity[2] > 1e4) {
    std::cerr << "Vellocity too high\n";
    exit(-41);
  }

  timer.registerTime("Gradient");
  timer.evaluateComponentsTime();
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
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        tempPhi[i][j][k] = _phi[i][j][k];
      }

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
        int cellSign = (cellPhi >= 0) ? 1 : -1;
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
          // cellSign *std::min(std::fabs(tempPhi[i][j][k]), theta *
          // _h.x());
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
  double maxDistance = 5 * _h[0];
  while (!cellsQueue.empty()) {
    int cellId = cellsQueue.top().second;
    cellsQueue.pop();
    int k = cellId / (_resolution.x() * _resolution.y());
    int j = (cellId % (_resolution.x() * _resolution.y())) / _resolution.x();
    int i = (cellId % (_resolution.x() * _resolution.y())) % _resolution.x();
    // TODO: check this caculation

    if (!processed[i][j][k]) {
      processed[i][j][k] = true;
      int distanceSignal = (_phi[i][j][k] >= 0) ? 1 : -1;
      if (std::fabs(_phi[i][j][k]) < 1e-7)
        distanceSignal = 1;
      double newPhi = _solveEikonal(glm::ivec3(i, j, k));

      if (newPhi < std::fabs(tempPhi[i][j][k]))
        tempPhi[i][j][k] = distanceSignal * newPhi;
    }
    if (tempPhi[i][j][k] > maxDistance)
      continue;
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
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[1]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[3]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[2]));
    mesh.addFace(face);
    break;
  case 14:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[1]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[2]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[0], vertices[3]));
    mesh.addFace(face);
    break;
  case 2:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[2]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[3]));
    mesh.addFace(face);
    break;
  case 13:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[3]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[1], vertices[2]));
    mesh.addFace(face);
    break;
  case 4:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[3]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[1]));
    mesh.addFace(face);
    break;
  case 11:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[1]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[2], vertices[3]));
    mesh.addFace(face);
    break;
  case 8:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[1]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[2]));
    mesh.addFace(face);
    break;
  case 7:
    face[0] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[0]));
    face[1] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[2]));
    face[2] = mesh.addVertex(_findSurfaceCoordinate(vertices[3], vertices[1]));
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
  if (std::isnan(phi1) || std::isnan(phi2) || std::isinf(phi1) ||
      std::isinf(phi2)) {
    std::stringstream message;
    message << "NOT VALID VALUE " << phi1 << " " << phi2 << "\n";
    throw(message.str().c_str());
  }
  return coord1 - phi1 * direction / (phi2 - phi1);
}

Eigen::Array3d LevelSet3::_findSurfaceCoordinate(Eigen::Array3i v1,
                                                 Eigen::Array3i v2) {
  Eigen::Array3d coord1, coord2, direction;
  Eigen::Array3d h(_h[0], _h[1], _h[2]);

  coord1 = v1.cast<double>().cwiseProduct(h) + h / 2.0;
  coord2 = v2.cast<double>().cwiseProduct(h) + h / 2.0;
  direction = coord2 - coord1;

  float phi1, phi2;
  phi1 = _phi[v1[0]][v1[1]][v1[2]];
  phi2 = _phi[v2[0]][v2[1]][v2[2]];

  if (std::isnan(phi1) || std::isnan(phi2) || std::isinf(phi1) ||
      std::isinf(phi2)) {
    std::stringstream message;
    message << "NOT VALID VALUE " << phi1 << " " << phi2 << "\n";
    throw(message.str().c_str());
  }

  return coord1 - phi1 * direction / (phi2 - phi1);
}

} // namespace Ramuh