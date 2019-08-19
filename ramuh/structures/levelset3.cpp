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
#include <blas/interpolator.h>
#include <blas/weno.h>

namespace Ramuh {

LevelSet3::LevelSet3() : LevelSet3(Eigen::Array3i(0, 0, 0)) {}

LevelSet3::LevelSet3(Eigen::Array3i resolution) : RegularGrid3(resolution) {
  _phi.resize(_resolution[0]);
  for (auto &row : _phi) {
    row.resize(_resolution[1]);
    for (auto &depth : row)
      depth.resize(_resolution[2], 1e8);
  }

  _isPressureSecondOrder = true;
}

void LevelSet3::setPressureSecondOrder(bool value) {
  _isPressureSecondOrder = value;
}

void LevelSet3::setResolution(Vector3i newResolution) {
  RegularGrid3::setResolution(newResolution);
  _phi.resize(_resolution[0]);
  for (auto &row : _phi) {
    row.resize(_resolution[1]);
    for (auto &depth : row)
      depth.resize(_resolution[2], 1e6);
  }
  _velocity.resize(_resolution[0] + 1);
  for (auto &row : _velocity) {
    row.resize(_resolution[1] + 1);
    for (auto &depth : row)
      depth.resize(_resolution[2] + 1);
  }
}

void LevelSet3::printVertexVelocity() {
  std::cerr << "==== x component: \n";
  for (int j = 0; j < _resolution[1] + 1; j++) {
    for (int i = 0; i < _resolution[0] + 1; i++) {
      std::cerr << _velocity[i][j][0][0] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "==== y component: \n";
  for (int j = 0; j < _resolution[1] + 1; j++) {
    for (int i = 0; i < _resolution[0] + 1; i++) {
      std::cerr << _velocity[i][j][0][1] << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
}

void LevelSet3::addSphereSurface(Vector3d center, double radius) {
  addSphereSurface(Eigen::Array3d(center[0], center[1], center[2]), radius);
}

void LevelSet3::addSphereSurface(Eigen::Array3d center, double radius) {
  for (int k = 0; k < _resolution[2]; k++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int i = 0; i < _resolution[0]; i++) {
        Eigen::Array3d h(_h[0], _h[1], _h[2]);
        Eigen::Array3d position =
            Eigen::Array3d(i, j, k).cwiseProduct(h) + h / 2.0;
        double distance = (position - center).matrix().norm() - radius;
        _phi[i][j][k] = std::min(_phi[i][j][k], distance);
      }
  checkCellMaterial();
}

void LevelSet3::addCubeSurface(Eigen::Array3d lower, Eigen::Array3d upper) {
  addCubeSurface(Vector3d(lower[0], lower[1], lower[2]),
                 Vector3d(upper[0], upper[1], upper[2]));
}

void LevelSet3::addCubeSurface(Vector3d lower, Vector3d upper) {
  Vector3d center = (lower + upper) / 2;
  Vector3d size = upper - lower;

  for (int k = 0; k < _resolution[2]; k++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int i = 0; i < _resolution[0]; i++) {
        Vector3d position(Vector3i(i, j, k) * _h + _h / 2);
        Vector3d distanceLower = (position - lower).abs();
        Vector3d distanceUpper = (position - upper).abs();
        double phi;
        if (position > lower && position < upper) {
          // Inside cube
          phi = std::min(std::fabs(_phi[i][j][k]),
                         std::min(distanceLower.min(), distanceUpper.min()));
          _phi[i][j][k] = -std::min(std::fabs(_phi[i][j][k]), phi);
        } else {
          // Compute distance to the cube
          position = position - center;
          position = position * 2 / size;

          position = position.abs();

          phi = 0.0;
          phi += std::max(0.0, position[0] - 1);
          phi += std::max(0.0, position[1] - 1);
          phi += std::max(0.0, position[2] - 1);
          phi = sqrt(phi);
          _phi[i][j][k] = std::min((_phi[i][j][k]), phi);
        }
      }
  checkCellMaterial();
}

void LevelSet3::checkCellMaterial() {
  _fluidCells.clear();
  _surfaceCells.clear();

// Mark all cells and faces as air
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];
    _material[i][j][k] = Material::FluidMaterial::AIR;
    _uFaceMaterial[i][j][k] = Material::FluidMaterial::AIR;
    _vFaceMaterial[i][j][k] = Material::FluidMaterial::AIR;
    _wFaceMaterial[i][j][k] = Material::FluidMaterial::AIR;
    _uFaceMaterial[i + 1][j][k] = Material::FluidMaterial::AIR;
    _vFaceMaterial[i][j + 1][k] = Material::FluidMaterial::AIR;
    _wFaceMaterial[i][j][k + 1] = Material::FluidMaterial::AIR;
  }

// Walk through all cell/faces again to mark them as fluid
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];
    _phi[i][j][k] = -1e-10;
    if (_phi[i][j][k] <= 0) {
      _uFaceMaterial[i + 1][j][k] = Material::FluidMaterial::FLUID;
      _vFaceMaterial[i][j + 1][k] = Material::FluidMaterial::FLUID;
      _wFaceMaterial[i][j][k + 1] = Material::FluidMaterial::FLUID;
#pragma omp critical
      {
        if (!_hasOppositeNeighborsWithMaterial(_id,
                                               Material::FluidMaterial::FLUID))
          _surfaceCells.emplace_back(_id);
        _fluidCells.emplace_back(_id);
      }
    }
  }
} // namespace Ramuh

void LevelSet3::addImplicitFunction() { NOT_IMPLEMENTED(); }

void LevelSet3::integrateLevelSet() {
  // std::vector<std::vector<double>> oldPhi;
  Matrix3<double> oldPhi;
  oldPhi.changeSize(_resolution);

#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    Eigen::Array3d position, backPosition, velocity, h(_h[0], _h[1], _h[2]),
        cellCenter;
    Eigen::Array3i index;
    double newPhi = _phi[i][j][k];

    position = Eigen::Array3d(i, j, k).cwiseProduct(h) +
               Eigen::Array3d(h[0] / 2, h[1] / 2, h[2] / 2);

    velocity[0] = _interpolateVelocityU(position);
    velocity[1] = _interpolateVelocityV(position);
    velocity[2] = _interpolateVelocityW(position);

    backPosition = position - (velocity * _dt);
    index = Eigen::floor(backPosition.cwiseQuotient(h)).cast<int>();

    // Check if inside domain
    if ((velocity.matrix().norm() > _tolerance) && index[0] >= 0 &&
        index[1] >= 0 && index[2] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1] && index[2] < _resolution[2]) {
      newPhi = _interpolatePhi(backPosition);
    }
    oldPhi[i][j][k] = newPhi;
  }
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];
    _phi[i][j][k] = oldPhi[i][j][k];
  }
}

void LevelSet3::macCormackAdvection() {
  Matrix3<double> newPhi;
  Matrix3<bool> phiChanged;
  newPhi.changeSize(_resolution);
  phiChanged.changeSize(_resolution, false);
  Eigen::Array3i resolution(_resolution[0], _resolution[1], _resolution[2]);
#pragma omp parallel for
  for (int it = 0; it < cellCount(); it++) {
    int _id = it;
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    if (std::fabs(_phi[i][j][k]) > 5 * _h[0]) {
      newPhi[i][j][k] = _phi[i][j][k];
      continue;
    }
    Eigen::Array3d position, h(_h[0], _h[1], _h[2]);
    Eigen::Vector3d velocity;
    Eigen::Array3i index;
    // Find threshold values for clamping
    double clamp[2]; // 0: min, 1: max
    clamp[0] = 1e8;
    clamp[1] = -1e8;

    // Find the point as if it was in previous timestep and work with it
    // Semi lagrangian advection
    position = Eigen::Array3d(i, j, k).cwiseProduct(h) + (h / 2.0);
    velocity[0] = _interpolateVelocityU(position);
    velocity[1] = _interpolateVelocityV(position);
    velocity[2] = _interpolateVelocityW(position);
    position -= velocity.array() * _dt;

    index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
    index[0] = std::max(0, std::min(resolution[0] - 1, index[0]));
    index[1] = std::max(0, std::min(resolution[0] - 1, index[1]));
    index[2] = std::max(0, std::min(resolution[0] - 1, index[2]));
    // Check if the index of the position is not a surface
    // bool isSurface = false;
    // int id = ijkToId(index[0], index[1], index[2]);
    // if (std::find(_surfaceCells.begin(), _surfaceCells.end(), _id) !=
    //         _surfaceCells.end() ||
    //     std::find(_surfaceCells.begin(), _surfaceCells.end(), id) !=
    //         _surfaceCells.end() ||
    //     _material[index[0]][index[1]][index[2]] ==
    //         Material::FluidMaterial::AIR) {
    //   isSurface = true;
    // }

    // TODO: look for a better solution
    if ((velocity.norm() > _tolerance) && index[0] >= 0 && index[1] >= 0 &&
        index[2] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1] && index[2] < _resolution[2]) {

      double phi_n;
      phi_n = _interpolatePhi(position, clamp[0], clamp[1]);

      // if (isSurface) {
      //   newPhi[i][j][k] = std::max(clamp[0], std::min(clamp[1], phi_n));
      // } else {
      // Advect forward in time the value at the position
      velocity[0] = _interpolateVelocityU(position);
      velocity[1] = _interpolateVelocityV(position);
      velocity[2] = _interpolateVelocityW(position);
      position += velocity.array() * _dt;

      try {
        double phi_n1_hat = _interpolatePhi(position);
        // Analyse the error from original to the forward advected
        double error = 0.5 * (_phi[i][j][k] - phi_n1_hat);
        newPhi[i][j][k] = std::max(clamp[0], std::min(clamp[1], phi_n + error));
      } catch (const char *error) {
        std::cerr << error;
        std::cerr << "Velocity " << velocity.transpose() << std::endl;
        std::cerr << "Original phi: " << _phi[i][j][k] << std::endl;
        std::cerr << "Interpolated phi_n: " << phi_n << std::endl;
        throw("levelSet3::macComarckAdvection\n");
      }
      if (!(newPhi[i][j][k] >= clamp[0] && newPhi[i][j][k] <= clamp[1]))
        newPhi[i][j][k] = phi_n;
    }
    // }
    // else if (isSurface)
    //   newPhi[i][j][k] = _interpolatePhi(position);
    else
      newPhi[i][j][k] = _phi[i][j][k];
  }
#pragma omp parallel for
  for (int i = 0; i < _resolution[0]; i++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int k = 0; k < _resolution[2]; k++)
        _phi[i][j][k] = newPhi[i][j][k];
}

void LevelSet3::advectWeno() {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Matrix3<double> newPhi;
  newPhi.changeSize(_resolution);

#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    std::vector<double> values(6);
    Eigen::Array3i ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];

    Eigen::Vector3d dPhi(0., 0., 0.);
    Eigen::Vector3d velocity;
    velocity[0] = _u[i][j][k];
    velocity[1] = _v[i][j][k];
    velocity[2] = _w[i][j][k];

    bool isNegative = true; // Velocity direction
    bool isShock = false;
    // X dimension
    if (velocity[0] >= 0)
      isNegative = false;
    if ((!isNegative && i > 0) || (isNegative && i == _resolution[0] - 1)) {
      isNegative = false;
      for (int ival = 0, ii = -3; ival < 7; ival++, ii++) {
        int index = std::min(_resolution[0] - 1, std::max(0, i + ii));
        double flux;
        flux = (_phi[index][j][k] +
                _phi[std::min(index + 1, _resolution[0] - 1)][j][k]) /
               2.0;
        flux *= _u[index][j][k];
        values[ival] = flux;
      }
    } else {
      isNegative = true;
      for (int ival = 0, ii = 3; ival < 7; ival++, ii--) {
        int index = std::min(_resolution[0] - 1, std::max(0, i + ii));
        double flux;
        flux = (_phi[index][j][k] + _phi[std::max(index - 1, 0)][j][k]) / 2.0;
        // flux *= _u[index][j][k];
        // flux = _phi[index][j][k];
        // flux *= (_u[index][j][k] +
        //          _u[std::min(index + 1, _resolution[0] - 1)][j][k]) /
        //         2.0;
        values[ival] = flux;
      }
    }
    dPhi[0] = Weno::evaluate(values, h[0], isNegative);
    dPhi[0] +=
        Weno::evaluate(std::vector<double>(values.begin() + 1, values.end()),
                       h[0], isNegative);
    // dPhi[0] /= 2.;

    // Y dimension
    isNegative = true;
    // isShock = false;
    if (velocity[1] >= 0)
      isNegative = false;
    // if (!isNegative) {
    if ((!isNegative && j > 0) || (isNegative && j == _resolution[1] - 1)) {
      isNegative = false;
      for (int ival = 0, jj = -3; ival < 7; ival++, jj++) {
        int index = std::min(_resolution[1] - 1, std::max(0, j + jj));
        double flux;
        flux = (_phi[i][index][k] +
                _phi[i][std::min(index + 1, _resolution[1] - 1)][k]) /
               2.0;
        flux *= _v[i][index][k];
        values[ival] = flux;
      }
    } else {
      isNegative = true;
      for (int ival = 0, jj = 3; ival < 7; ival++, jj--) {
        int index = std::min(_resolution[1] - 1, std::max(0, j + jj));
        double flux;
        flux = (_phi[i][index][k] +
                _phi[i][std::min(index + 1, _resolution[1] - 1)][k]) /
               2.0;
        flux *= _v[i][index][k];
        values[ival] = flux;
      }
    }
    dPhi[1] = 0; // Weno::evaluate(values, h[0], isNegative);

    // Z dimension
    isNegative = true;
    // isShock = false;
    if (velocity[2] >= 0)
      isNegative = false;
    // if (!isNegative) {
    if ((!isNegative && k > 0) || (isNegative && k == _resolution[2] - 1)) {
      isNegative = false;
      for (int ival = 0, kk = -3; ival < 7; ival++, kk++) {
        int index = std::min(_resolution[2] - 1, std::max(0, k + kk));
        double flux;
        flux = (_phi[i][j][index] +
                _phi[i][j][std::min(index + 1, _resolution[2] - 1)]) /
               2.0;
        flux *= _w[i][j][index];
        values[ival] = flux;
      }
    } else {
      isNegative = true;
      for (int ival = 0, kk = 3; ival < 7; ival++, kk--) {
        int index = std::min(_resolution[2] - 1, std::max(0, k + kk));
        double flux;
        flux = (_phi[i][j][index] +
                _phi[i][j][std::min(index + 1, _resolution[2] - 1)]) /
               2.0;
        flux *= _w[i][j][index];
        values[ival] = flux;
      }
    }
    dPhi[2] = 0; // Weno::evaluate(values, h[0], isNegative);

    // Proceed to time integration and material derivative
    // Euler method
    if (i > 0 && _phi[i - 1][j][k] < 0)
      int a = 0;
    newPhi[i][j][k] = -dPhi.dot(velocity) * _dt + _phi[i][j][k];
  }

#pragma omp parallel for
  for (size_t id = 0; id < cellCount(); id++) {
    Eigen::Array3i ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    _phi[i][j][k] = newPhi[i][j][k];
  }
}

double LevelSet3::_interpolatePhi(Eigen::Array3d position) {
  double _min, _max;
  return _interpolatePhi(position, _min, _max);
}

double LevelSet3::_interpolatePhi(Eigen::Array3d position, double &_min,
                                  double &_max) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution[0], _resolution[1], _resolution[2]);
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>().cwiseProduct(h) + h / 2.0;

  if (position[0] < cellCenter[0])
    index[0]--;
  if (position[1] < cellCenter[1])
    index[1]--;
  if (position[2] < cellCenter[2])
    index[2]--;

  std::vector<int> iCandidates, jCandidates, kCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
    kCandidates.push_back(index[2] + i);
  }

  std::vector<Eigen::Array3d> points;
  std::vector<double> values;

  _min = 1e8;
  _max = -1e8;

  for (auto w : kCandidates)
    for (auto v : jCandidates)
      for (auto u : iCandidates) {
        int x, y, z;
        points.emplace_back((u + 0.5) * h[0], (v + 0.5) * h[1],
                            (w + 0.5) * h[2]);
        x = std::max(0, std::min(resolution[0] - 1, u));
        y = std::max(0, std::min(resolution[1] - 1, v));
        z = std::max(0, std::min(resolution[2] - 1, w));
        values.emplace_back(_phi[x][y][z]);

        _min = std::min(values.back(), _min);
        _max = std::max(values.back(), _max);
      }
  double phi = Interpolator::tricubic(position, points, values);
  return std::min(_max, std::max(_min, phi));
  // return phi;
}

double LevelSet3::_interpolatePhi(Eigen::Array3d position, int originalSignal) {
  Eigen::Array3d h(_h[0], _h[1], _h[2]);
  Eigen::Array3i resolution(_resolution[0], _resolution[1], _resolution[2]);
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>().cwiseProduct(h) + h / 2.0;
  index[0] = std::max(0, std::min(_resolution[0], index[0]));
  index[1] = std::max(0, std::min(_resolution[1], index[1]));
  index[2] = std::max(0, std::min(_resolution[2], index[2]));
  std::vector<int> iCandidates, jCandidates, kCandidates;

  if (index[0] >= 0 && index[0] < _resolution[0]) {
    iCandidates.push_back(index[0]);
  }
  if (index[1] >= 0 && index[1] < _resolution[1]) {
    jCandidates.push_back(index[1]);
  }
  if (index[2] >= 0 && index[2] < _resolution[2]) {
    kCandidates.push_back(index[2]);
  }
  if (position[0] > cellCenter[0] && index[0] < resolution[0] - 1) {
    iCandidates.push_back(index[0] + 1);
  } else if (position[0] < cellCenter[0] && index[0] > 0) {
    iCandidates.push_back(index[0] - 1);
  }
  if (position[1] > cellCenter[1] && index[1] < resolution[1] - 1) {
    jCandidates.push_back(index[1] + 1);
  } else if (position[1] < cellCenter[1] && index[1] > 0) {
    jCandidates.push_back(index[1] - 1);
  }
  if (position[2] > cellCenter[2] && index[2] < resolution[2] - 1) {
    kCandidates.push_back(index[2] + 1);
  } else if (position[2] < cellCenter[2] && index[2] > 0) {
    kCandidates.push_back(index[2] - 1);
  }

  // Catmull-Rom like interpolation
  double newPhi = 0., distance = 0., distanceCount = 0.;
  Eigen::Array2d clamp; // 0: min, 1: max
  clamp[0] = 1e8;
  clamp[1] = -1e8;
  for (auto u : iCandidates)
    for (auto v : jCandidates)
      for (auto w : kCandidates) {
        if (u < 0 || u >= _resolution[0] || v < 0 || v >= _resolution[1] ||
            w < 0 || w >= _resolution[2])
          continue;
        Eigen::Array3d centerPosition = Eigen::Array3d(u, v, w) * h + h / 2.0;
        distance = (position - centerPosition).matrix().norm();
        if (distance < 1e-6)
          return _phi[u][v][w];
        double phiValue = _phi[u][v][w];
        int phiSignal = (phiValue >= 0) ? 1 : -1;
        if (phiSignal == originalSignal) {
          distanceCount += 1. / distance;
          newPhi += _phi[u][v][w] / distance;
          clamp[0] = std::min(clamp[0], _phi[u][v][w]);
          clamp[1] = std::max(clamp[1], _phi[u][v][w]);
        }
      }
  if (distanceCount == 0 || distanceCount > 1e8) {
    std::cerr << "LevelSet3::interpolatePhi: distanceCount error: "
              << distanceCount << "\n";
    std::cerr << "center: " << cellCenter.transpose() << std::endl;
    std::cerr << "position: " << position.transpose() << std::endl;
    std::cerr << "phiSignal: " << originalSignal << std::endl;
    throw("LevelSet3::interpolatePhi: distanceCount zero or infinity\n");
  }
  // return std::max(clamp[0], std::min(clamp[1], newPhi / distanceCount));
  return newPhi / distanceCount;
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

  // Solve pressure Poisson equation
  divergent = Eigen::VectorXd::Zero(cellCount());
// int ghostId = _getMapId(-1);
#pragma omp parallel for
  for (int it = 0; it < _fluidCells.size(); it++) {
    int _id = _fluidCells[it];
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    if (_phi[i][j][k] > 0.0)
      continue;
    int id = _getMapId(_id);

    int validCells = 0;
    double h2 = _h[0] * _h[0];
    double centerWeight = 0.0;
    std::vector<Eigen::Triplet<double>> threadTriplet;
    threadTriplet.clear();
    Eigen::Array3d center = Eigen::Array3d(i, j, k).cwiseProduct(h) + h / 2.0;

    if (i > 0) {
      validCells++;
      if (_phi[i - 1][j][k] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i - 1, j, k)),
                                   -1 / (h2));
      } else {
        if (_isPressureSecondOrder) {
          double theta;
          // Eigen::Array3d surface = _findSurfaceCoordinate(
          // Eigen::Array3i(i, j, k), Eigen::Array3i(i - 1, j, k));
          // theta = (surface - center).matrix().norm() / _h[0];
          theta = -_phi[i][j][k] / (_phi[i - 1][j][k] - _phi[i][j][k]);
          centerWeight += (1 / (h2 * theta));
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (i < _resolution[0] - 1) {
      validCells++;
      if (_phi[i + 1][j][k] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i + 1, j, k)),
                                   -1 / (h2));
      } else {
        if (_isPressureSecondOrder) {
          double theta;
          // Eigen::Array3d surface = _findSurfaceCoordinate(
          // Eigen::Array3i(i, j, k), Eigen::Array3i(i + 1, j, k));
          // theta = (surface - center).matrix().norm() / _h[0];
          theta = -_phi[i][j][k] / (_phi[i + 1][j][k] - _phi[i][j][k]);
          centerWeight += (1 / (h2 * theta));
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (j > 0) {
      validCells++;
      if (_phi[i][j - 1][k] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j - 1, k)),
                                   -1 / (h2));
      } else {
        if (_isPressureSecondOrder) {
          double theta;
          // Eigen::Array3d surface = _findSurfaceCoordinate(
          // Eigen::Array3i(i, j, k), Eigen::Array3i(i, j - 1, k));
          // theta = (surface - center).matrix().norm() / _h[1];
          theta = -_phi[i][j][k] / (_phi[i][j - 1][k] - _phi[i][j][k]);
          centerWeight += (1 / (h2 * theta));
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (j < _resolution[1] - 1) {
      validCells++;
      if (_phi[i][j + 1][k] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j + 1, k)),
                                   -1 / (h2));
      } else {
        if (_isPressureSecondOrder) {
          double theta;
          // Eigen::Array3d surface = _findSurfaceCoordinate(
          // Eigen::Array3i(i, j, k), Eigen::Array3i(i, j + 1, k));
          // theta = (surface - center).matrix().norm() / _h[1];
          theta = -_phi[i][j][k] / (_phi[i][j + 1][k] - _phi[i][j][k]);
          centerWeight += (1 / (h2 * theta));
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (k > 0) {
      validCells++;
      if (_phi[i][j][k - 1] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j, k - 1)),
                                   -1 / (h2));
      } else {
        if (_isPressureSecondOrder) {
          double theta;
          // Eigen::Array3d surface = _findSurfaceCoordinate(
          // Eigen::Array3i(i, j, k), Eigen::Array3i(i, j, k - 1));
          // theta = (surface - center).matrix().norm() / _h[2];
          theta = -_phi[i][j][k] / (_phi[i][j][k - 1] - _phi[i][j][k]);
          centerWeight += (1 / (h2 * theta));
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    if (k < _resolution[2] - 1) {
      validCells++;
      if (_phi[i][j][k + 1] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijkToId(i, j, k + 1)),
                                   -1 / (h2));
      } else {
        if (_isPressureSecondOrder) {
          double theta;
          // Eigen::Array3d surface = _findSurfaceCoordinate(
          // Eigen::Array3i(i, j, k), Eigen::Array3i(i, j, k + 1));
          // theta = (surface - center).matrix().norm() / _h[2];
          theta = -_phi[i][j][k] / (_phi[i][j][k + 1] - _phi[i][j][k]);
          centerWeight += (1 / (h2 * theta));
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          // threadTriplet.emplace_back(id, ghostId, -1 / (h2));
        }
      }
    }
    threadTriplet.emplace_back(id, id, centerWeight);
#pragma omp critical
    {
      triplets.insert(triplets.end(), threadTriplet.begin(),
                      threadTriplet.end());
    }

    divergent[id] = 0;
    divergent[id] -= (_u[i + 1][j][k] - _u[i][j][k]) / _h[0];
    divergent[id] -= (_v[i][j + 1][k] - _v[i][j][k]) / _h[1];
    divergent[id] -= (_w[i][j][k + 1] - _w[i][j][k]) / _h[2];
    divergent[id] /= _dt;
  }
  // divergent[ghostId] = 0;
  // triplets.emplace_back(ghostId, ghostId, 1);

  int nCells = idMap.size();
  pressure = Eigen::VectorXd::Zero(nCells);
  divergent.conservativeResize(nCells);

  Eigen::SparseMatrix<double> pressureMatrix(nCells, nCells);
  pressureMatrix.setFromTriplets(triplets.begin(), triplets.end());

  // SOlve pressure Poisson system
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.setTolerance(1e-10);
  // solver.setMaxIterations(100);
  solver.compute(pressureMatrix);
  pressure = solver.solve(divergent);
  // std::cerr << divergent.transpose() << std::endl;

  // std::cerr << "Solved with " << solver.iterations() << " iterations.\n";
  // std::cerr << "Solver error: " << solver.error() << std::endl;
  // {
  //   Eigen::VectorXd x = pressureMatrix * pressure;
  //   std::cerr << "L2 norm: " << (x - divergent).norm() << std::endl;
  //   if (!x.isApprox(divergent, 1e-2))
  //     std::cerr << "Error is too high\n";
  // }

  // Correct velocity through pressure gradient
  _maxVelocity[0] = _maxVelocity[1] = _maxVelocity[2] = -1e8;
#pragma omp parallel for
  // TODO: Change to fluid FACES instead of CELLS
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
    double weight;
    auto it = idMap.find(id);

    pressure1 = (it == idMap.end()) ? 0.0 : pressure[it->second];
    pressure2 = 0.0;
    if (i > 0) {
      int neighId = ijkToId(i - 1, j, k);
      auto it = idMap.find(ijkToId(i - 1, j, k));
      if (it == idMap.end())
        pressure2 = 0.0;
      else
        pressure2 = pressure[it->second];
      if (_isPressureSecondOrder &&
          std::signbit(_phi[i][j][k]) != std::signbit(_phi[i - 1][j][k])) {
        if (_phi[i][j][k] <= 0)
          weight =
              -_phi[i][j][k] / std::fabs(_phi[i - 1][j][k] - _phi[i][j][k]);
        else
          weight =
              -_phi[i - 1][j][k] / std::fabs(_phi[i - 1][j][k] - _phi[i][j][k]);
        if (std::fabs(weight) < 1e-10) {
          weight = 1;
          pressure1 = 0;
        }
      } else
        weight = 1.0;
      _u[i][j][k] -= (_dt * (pressure1 - pressure2) / (_h[0] * weight));
      if (std::isnan(_u[i][j][k]) || std::isinf(_u[i][j][k])) {
        std::cerr << "computePressure: Infinite u velocity component\n";
        std::cerr << "Pressures: " << pressure1 << " - " << pressure2
                  << std::endl;
        std::cerr << "Weight: " << weight << ": " << _phi[i][j][k] << std::endl;
      }
      _maxVelocity[0] = std::max(_maxVelocity[0], std::fabs(_u[i][j][k]));
    }
    pressure2 = 0.0;
    int neighId = ijkToId(i, j - 1, k);
    if (j > 0) {
      auto it = idMap.find(ijkToId(i, j - 1, k));
      if (it == idMap.end())
        pressure2 = 0.0;
      else
        pressure2 = pressure[it->second];
      if (_isPressureSecondOrder &&
          std::signbit(_phi[i][j][k]) != std::signbit(_phi[i][j - 1][k])) {
        if (_phi[i][j][k] <= 0)
          weight =
              -_phi[i][j][k] / std::fabs(_phi[i][j - 1][k] - _phi[i][j][k]);
        else
          weight =
              -_phi[i][j - 1][k] / std::fabs(_phi[i][j - 1][k] - _phi[i][j][k]);
        if (std::fabs(weight) < 1e-10) {
          weight = 1;
          pressure1 = 0;
        }
      } else
        weight = 1.0;
      _v[i][j][k] -= (_dt * (pressure1 - pressure2) / (_h[1] * weight));
      if (std::isnan(_v[i][j][k]) || std::isinf(_v[i][j][k])) {
        std::cerr << "computePressure: Infinite v velocity component";
        std::cerr << "Pressures: " << pressure1 << " - " << pressure2
                  << std::endl;
        std::cerr << "Weight: " << weight << ": " << _phi[i][j][k] << std::endl;
      }
      _maxVelocity[1] = std::max(_maxVelocity[1], std::fabs(_v[i][j][k]));
    }
    pressure2 = 0.0;
    if (k > 0) {
      auto it = idMap.find(ijkToId(i, j, k - 1));
      if (it == idMap.end())
        pressure2 = 0.0;
      else
        pressure2 = pressure[it->second];
      if (_isPressureSecondOrder &&
          std::signbit(_phi[i][j][k]) != std::signbit(_phi[i][j][k - 1])) {
        if (_phi[i][j][k] <= 0)
          weight =
              -_phi[i][j][k] / std::fabs(_phi[i][j][k - 1] - _phi[i][j][k]);
        else
          weight =
              -_phi[i][j][k - 1] / std::fabs(_phi[i][j][k - 1] - _phi[i][j][k]);
        if (std::fabs(weight) < 1e-10) {
          weight = 1;
          pressure1 = 0;
        }
      } else
        weight = 1.0;
      _w[i][j][k] -= (_dt * (pressure1 - pressure2) / (_h[2] * weight));
      if (std::isnan(_w[i][j][k]) || std::isinf(_w[i][j][k])) {
        std::cerr << "computePressure: Infinite w velocity component";
        std::cerr << "Pressures: " << pressure1 << " - " << pressure2
                  << std::endl;
        std::cerr << "Weight: " << weight << ": " << _phi[i][j][k] << std::endl;
      }
      _maxVelocity[2] = std::max(_maxVelocity[2], std::fabs(_w[i][j][k]));
    }
  }
  if (_maxVelocity[0] > 1e4 || _maxVelocity[1] > 1e4 || _maxVelocity[2] > 1e4) {
    std::cerr << "Vellocity too high: " << _maxVelocity[0] << _maxVelocity[1]
              << _maxVelocity[2] << std::endl;
    exit(-41);
  }

  {
    std::ofstream file;
    file.open("results/pressValues");
    if (file.is_open()) {
      for (int id = 0; id < cellCount(); id++) {
        auto it = idMap.find(id);
        if (it == idMap.end()) {
          file << " 0 ";
        } else {
          file << " " << pressure[it->second] << " ";
        }
      }
    }
    std::cout << "pressure written\n";
    file.close();
    file.open("results/divergent");

    for (int id = 0; id < cellCount(); id++) {
      auto it = idMap.find(id);
      if (it == idMap.end()) {
        file << " 0 ";
      } else {
        file << " " << divergent[it->second] << " ";
      }
    }
    file.close();
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
               std::fabs(_phi[std::min(_resolution[0] - 1, i + 1)][j][k]));
  distances[1] =
      std::min(std::fabs(_phi[i][std::max(0, j - 1)][k]),
               std::fabs(_phi[i][std::min(_resolution[1] - 1, j + 1)][k]));
  distances[2] =
      std::min(std::fabs(_phi[i][j][std::max(0, k - 1)]),
               std::fabs(_phi[i][j][std::min(_resolution[2] - 1, k + 1)]));

  std::sort(distances.begin(), distances.end());

  double newPhi = distances[0] + _h[0];
  if (newPhi > distances[1]) {
    double h2 = _h[0] * _h[0];
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
  if (std::fabs(newPhi) < _tolerance || std::isnan(newPhi) || newPhi == 0 ||
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

  // Copying a backup values from previous frmae
  tempPhi.changeSize(_resolution, 1e8);
  for (int i = 0; i < _resolution[0]; i++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int k = 0; k < _resolution[2]; k++) {
        tempPhi[i][j][k] = _phi[i][j][k];
      }

  processed.changeSize(_resolution, false);

#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array3i ijk = idToijk(_id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    glm::vec3 position = glm::vec3(i * _h[0], j * _h[1], k * _h[2]);
    glm::vec3 intersections[3];
    int nintersecs = 0;
    double cellPhi = _phi[i][j][k];
    int cellSign = (cellPhi >= 0) ? 1 : -1;
    auto signBit = std::signbit(cellSign);
    if (cellPhi == 0)
      cellSign = 1;

    bool isSurface = false, intersected = false;
    double theta = 1e8;
    // Find whether edge has an intersection with fluid surface
    if (i < _resolution[0] - 1 && signBit != std::signbit(_phi[i + 1][j][k])) {
      isSurface = true;
      intersected = true;
      theta =
          std::min(theta, std::fabs(cellPhi / (cellPhi - _phi[i + 1][j][k])));
      intersections[nintersecs] =
          glm::vec3(position[0] + theta * _h[0], position[1], position[2]);
    }
    if (i > 0 && signBit != std::signbit(_phi[i - 1][j][k])) {
      isSurface = true;
      intersected = true;
      theta =
          std::min(theta, std::fabs(cellPhi / (cellPhi - _phi[i - 1][j][k])));
      intersections[nintersecs] =
          glm::vec3(position[0] - theta * _h[0], position[1], position[2]);
    }
    theta = 1e8;
    if (intersected) {
      intersected = false;
      nintersecs++;
    }
    if (j < _resolution[1] - 1 && signBit != std::signbit(_phi[i][j + 1][k])) {
      isSurface = true;
      intersected = true;
      theta =
          std::min(theta, std::fabs(cellPhi / (cellPhi - _phi[i][j + 1][k])));
      intersections[nintersecs] =
          glm::vec3(position[0], position[1] + theta * _h[1], position[2]);
    }
    if (j > 0 && signBit != std::signbit(_phi[i][j - 1][k])) {
      isSurface = true;
      intersected = true;
      theta =
          std::min(theta, std::fabs(cellPhi / (cellPhi - _phi[i][j - 1][k])));
      intersections[nintersecs] =
          glm::vec3(position[0], position[1] - theta * _h[1], position[2]);
    }
    theta = 1e8;
    if (intersected) {
      intersected = false;
      nintersecs++;
    }
    if (k < _resolution[2] - 1 && signBit != std::signbit(_phi[i][j][k + 1])) {
      isSurface = true;
      intersected = true;
      theta =
          std::min(theta, std::fabs(cellPhi / (cellPhi - _phi[i][j][k + 1])));
      intersections[nintersecs] =
          glm::vec3(position[0], position[1], position[2] + theta * _h[2]);
    }
    if (k > 0 && signBit != std::signbit(_phi[i][j][k - 1])) {
      isSurface = true;
      intersected = true;
      theta =
          std::min(theta, std::fabs(cellPhi / (cellPhi - _phi[i][j][k - 1])));
      intersections[nintersecs] =
          glm::vec3(position[0], position[1], position[2] - theta * _h[2]);
    }
    if (intersected) {
      intersected = false;
      nintersecs++;
    }

    if (isSurface) {
      glm::vec3 proj;
      // Find intersection point (nearest to the intersecting geometry)
      if (nintersecs == 1) {
        proj = intersections[0];
      } else if (nintersecs == 2) {
        proj = Geometry::closestPointPlane(position, intersections[0],
                                           intersections[1]);
      } else if (nintersecs == 3) {
        proj = Geometry::closestPointTriangle(
            position, intersections[0], intersections[1], intersections[2]);
      }
      // Distance itself
      double distance = glm::length(proj - position);
      if (std::isnan(distance) || distance == 0 || distance > 5)
        float d = distance;
      if (std::fabs(distance) < _tolerance)
        distance = 0;

      // Keep signalfrom original distance
      tempPhi[i][j][k] = cellSign * distance;
      // _phi[i][j][k];
      // cellSign * _solveEikonal(glm::ivec3(i, j, k));
      // cellSign *std::min(std::fabs(tempPhi[i][j][k]), theta *
      // _h[0]);
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
    int k = cellId / (_resolution[0] * _resolution[1]);
    int j = (cellId % (_resolution[0] * _resolution[1])) / _resolution[0];
    int i = (cellId % (_resolution[0] * _resolution[1])) % _resolution[0];
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
    if (i < _resolution[0] - 1 &&
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
    if (j < _resolution[1] - 1 &&
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
    if (k < _resolution[2] - 1 &&
        cellsAdded.find(ijkToId(i, j, k + 1)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i][j][k + 1]), ijkToId(i, j, k + 1)));
      cellsAdded.insert(ijkToId(i, j, k + 1));
    }
  }

  for (int i = 0; i < _resolution[0]; i++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int k = 0; k < _resolution[2]; k++)
        _phi[i][j][k] = tempPhi[i][j][k];
}

void LevelSet3::writeLevelSetValue(const char *filename) {
  std::ofstream file;
  file.open(filename, std::ofstream::out | std::ofstream::trunc);

  if (!file.is_open()) {
    std::cerr << "Levelset3::writeLevelSetValue: Failed to open file "
              << filename << std::endl;
    return;
  }

  for (int i = 0; i < _resolution[0]; i++) {
    for (int j = 0; j < _resolution[1]; j++) {
      for (int k = 0; k < _resolution[2]; k++) {
        file << _phi[i][j][k] << ' ';
      }
      file << std::endl;
    }
    file << std::endl;
  }
}

void LevelSet3::readLevelSetValue(const char *filename) {
  std::ifstream file;
  file.open(filename, std::ifstream::in);

  if (!file.is_open()) {
    std::cerr << "Levelset3::readLevelSetValue: Failed to open file "
              << filename << std::endl;
    return;
  }
  for (int i = 0; i < _resolution[0]; i++) {
    for (int j = 0; j < _resolution[1]; j++) {
      for (int k = 0; k < _resolution[2]; k++) {
        file >> _phi[i][j][k];
      }
    }
  }
}

TriangleMesh LevelSet3::marchingTetrahedra() {
  TriangleMesh mesh;

  // Look for octants that have different levelset signal
  for (int i = -1; i < _resolution[0]; i++)
    for (int j = -1; j < _resolution[1]; j++)
      for (int k = -1; k < _resolution[2]; k++) {
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

  if (mesh.getVerticesSize() == 0) {
    std::cerr << "All fluid is gone!\n";
    throw("Marching tetrahedra: No fluid");
  }
  return mesh;
}

void LevelSet3::_triangulate(std::vector<glm::ivec3> vertices,
                             TriangleMesh &mesh) {
  // Check if vertices are in proper order and invert them if necessary
  glm::vec3 coords[4];
  double values[4];
  glm::vec3 h(_h[0], _h[1], _h[2]);
  for (int vert = 0; vert < 4; vert++) {
    coords[vert] = glm::vec3(vertices[vert]) * h + h / 2.0f;
    if ((vertices[vert][0] >= 0 && vertices[vert][0] < _resolution[0]) &&
        (vertices[vert][1] >= 0 && vertices[vert][1] < _resolution[1]) &&
        (vertices[vert][2] >= 0 && vertices[vert][2] < _resolution[2]))
      values[vert] =
          _phi[vertices[vert][0]][vertices[vert][1]][vertices[vert][2]];
    else {
      int i, j, k;
      i = std::max(0, std::min(_resolution[0] - 1, vertices[vert][0]));
      j = std::max(0, std::min(_resolution[1] - 1, vertices[vert][1]));
      k = std::max(0, std::min(_resolution[2] - 1, vertices[vert][2]));

      if (vertices[vert][0] < 0) {
        coords[vert][0] += h[0] / 2.f;
        values[vert] = (_phi[i][j][k] <= 0) ? -_phi[i][j][k] : _phi[i][j][k];
      }
      if (vertices[vert][1] < 0) {
        coords[vert][1] += h[1] / 2.f;
        values[vert] = (_phi[i][j][k] <= 0) ? -_phi[i][j][k] : _phi[i][j][k];
      }
      if (vertices[vert][2] < 0) {
        coords[vert][2] += h[2] / 2.f;
        values[vert] = (_phi[i][j][k] <= 0) ? -_phi[i][j][k] : _phi[i][j][k];
      }
      if (vertices[vert][0] > _resolution[0] - 1) {
        coords[vert][0] -= h[0] / 2.f;
        values[vert] = (_phi[i][j][k] <= 0) ? -_phi[i][j][k] : _phi[i][j][k];
      }
      if (vertices[vert][1] > _resolution[1] - 1) {
        coords[vert][1] -= h[1] / 2.f;
        values[vert] = (_phi[i][j][k] <= 0) ? -_phi[i][j][k] : _phi[i][j][k];
      }
      if (vertices[vert][2] > _resolution[2] - 1) {
        coords[vert][2] -= h[2] / 2.f;
        values[vert] = (_phi[i][j][k] <= 0) ? -_phi[i][j][k] : _phi[i][j][k];
      }
    }
  }

  if (glm::dot(glm::cross(coords[1] - coords[0], coords[2] - coords[0]),
               coords[3] - coords[0]) > 0) {
    auto aux = coords[1];
    coords[1] = coords[2];
    coords[2] = aux;

    auto aux2 = values[1];
    values[1] = values[2];
    values[2] = aux2;
  }

  _triangulate(coords, values, mesh);
}

void LevelSet3::_triangulate(glm::vec3 coords[4], double values[4],
                             TriangleMesh &mesh) {

  int triIndex = 0;
  if (values[0] < 0)
    triIndex |= 1;
  if (values[1] < 0)
    triIndex |= 2;
  if (values[2] < 0)
    triIndex |= 4;
  if (values[3] < 0)
    triIndex |= 8;

  glm::vec3 vertex;
  glm::ivec3 face;
  // 15 total cases to check
  switch (triIndex) {
  case 0:
  case 15:
    break;
  case 1:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[1], values[0], values[1]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[3], values[0], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[2], values[0], values[2]));
    mesh.addFace(face);
    break;
  case 14:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[1], values[0], values[1]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[2], values[0], values[2]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[3], values[0], values[3]));
    mesh.addFace(face);
    break;
  case 2:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[0], values[1], values[0]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[2], values[1], values[2]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[3], values[1], values[3]));
    mesh.addFace(face);
    break;
  case 13:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[0], values[1], values[0]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[3], values[1], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[2], values[1], values[2]));
    mesh.addFace(face);
    break;
  case 4:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[0], values[2], values[0]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[1], values[2], values[1]));
    mesh.addFace(face);
    break;
  case 11:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[0], values[2], values[0]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[1], values[2], values[1]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    mesh.addFace(face);
    break;
  case 8:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[3], coords[0], values[3], values[0]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[3], coords[1], values[3], values[1]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[3], coords[2], values[3], values[2]));
    mesh.addFace(face);
    break;
  case 7:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[3], coords[0], values[3], values[0]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[3], coords[2], values[3], values[2]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[3], coords[1], values[3], values[1]));
    mesh.addFace(face);
    break;
  case 3:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[2], values[0], values[2]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[3], values[1], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[3], values[0], values[3]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[2], values[1], values[2]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[3], values[1], values[3]));
    mesh.addFace(face);
    break;
  case 12:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[2], values[0], values[2]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[3], values[0], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[3], values[1], values[3]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[3], values[1], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[2], values[1], values[2]));
    mesh.addFace(face);
    break;
  case 5:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[1], values[0], values[1]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[3], values[0], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[2], values[1], values[2]));
    mesh.addFace(face);
    break;
  case 10:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[1], values[0], values[1]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[3], values[0], values[3]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[2], values[1], values[2]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    mesh.addFace(face);
    break;
  case 9:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[1], values[0], values[1]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[2], values[0], values[2]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[3], values[1], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    mesh.addFace(face);
    break;
  case 6:
    face[0] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[1], values[0], values[1]));
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[0], coords[2], values[0], values[2]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    mesh.addFace(face);
    face[1] = mesh.addVertex(
        _findSurfaceCoordinate(coords[2], coords[3], values[2], values[3]));
    face[2] = mesh.addVertex(
        _findSurfaceCoordinate(coords[1], coords[3], values[1], values[3]));
    mesh.addFace(face);
    break;
  }
}

// TODO: remove this function
glm::vec3 LevelSet3::_findSurfaceCoordinate(glm::vec3 coord1, glm::vec3 coord2,
                                            float phi1, float phi2) {
  glm::vec3 direction;
  direction = coord2 - coord1;

  if (phi1 == 0.0)
    return coord1;
  if (phi2 == 0.0)
    return coord2;

  if (std::signbit(phi1) == std::signbit(phi2)) {
    std::stringstream message;
    message << "LevelSet3::_findSurfaceCoordinate: There is no surface between "
               "points\n";
    throw(message.str().c_str());
  }
  if (std::isnan(phi1) || std::isnan(phi2) || std::isinf(phi1) ||
      std::isinf(phi2)) {
    std::stringstream message;
    message << "NOT VALID VALUE " << phi1 << " " << phi2 << "\n";
    throw(message.str().c_str());
  }
  return coord1 - (phi1 * direction) / (phi2 - phi1);
}

glm::vec3 LevelSet3::_findSurfaceCoordinate(glm::ivec3 v1, glm::ivec3 v2) {
  glm::vec3 coord1, coord2, direction;
  coord1 = glm::vec3(v1[0] * _h[0], v1[1] * _h[1], v1[2] * _h[2]);
  coord2 = glm::vec3(v2[0] * _h[0], v2[1] * _h[1], v2[2] * _h[2]);
  direction = coord2 - coord1;

  float phi1, phi2;
  phi1 = _phi[v1[0]][v1[1]][v1[2]];
  phi2 = _phi[v2[0]][v2[1]][v2[2]];
  if (phi1 == 0.0)
    return coord1;
  if (phi2 == 0.0)
    return coord2;

  if (std::signbit(phi1) == std::signbit(phi2)) {
    std::stringstream message;
    message << "LevelSet3::_findSurfaceCoordinate: There is no surface between "
               "points\n";
    throw(message.str().c_str());
  }
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

  if (phi1 == 0.0)
    return coord1;
  if (phi2 == 0.0)
    return coord2;

  if (std::signbit(phi1) == std::signbit(phi2)) {
    std::stringstream message;
    message << "LevelSet3::_findSurfaceCoordinate: There is no surface between "
               "points\n";
    throw(message.str().c_str());
  }

  if (std::isnan(phi1) || std::isnan(phi2) || std::isinf(phi1) ||
      std::isinf(phi2)) {
    std::stringstream message;
    message << "NOT VALID VALUE " << phi1 << " " << phi2 << "\n";
    throw(message.str().c_str());
  }

  return coord1 - phi1 * direction / (phi2 - phi1);
}

} // namespace Ramuh