#include <fluids/levelset_fluid3.h>
#include <blas/interpolator.h>
#include <blas/weno.h>

namespace Leviathan {

LevelSetFluid3::LevelSetFluid3()
    : LevelSetFluid3(Eigen::Array3i(32, 32, 32),
                     Ramuh::BoundingBox3(Eigen::Array3d(-1, -1, -1),
                                         Eigen::Array3d(1, 1, 1))) {}

LevelSetFluid3::LevelSetFluid3(Eigen::Array3i gridSize,
                               Ramuh::BoundingBox3 domain)
    : Ramuh::MacGrid3(domain, gridSize) {
  _phiId = newScalarLabel("phi", 1e8);
  _velocityId = newFaceScalarLabel("velocity");
  //  --> Even though velocities are vector, MAC grid split them into scalar
  //  pieces for each coordinate
  _dt = 1. / 60;
  _tolerance = 1e-10;

  newScalarLabel("newPhi");
}

void LevelSetFluid3::advectSemiLagrangean() {
  size_t nCells = cellCount();
  auto &newPhi = getScalarData("newPhi");
  auto &phi = getScalarData(_phiId);

#pragma omp parallel for
  for (size_t id = 0; id < nCells; id++) {
    Eigen::Array3d position, h, velocity;
    h = getH();

    velocity[0] = __interpolateVelocityU(position);
    velocity[1] = __interpolateVelocityV(position);
    velocity[2] = __interpolateVelocityW(position);

    position = position - (velocity * _dt);
    Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
    // Check if inside domain
    double newValue = phi[id];
    if ((velocity.matrix().norm() > _tolerance) && index[0] >= 0 &&
        index[1] >= 0 && index[2] >= 0 && index[0] < _gridSize[0] &&
        index[1] < _gridSize[1] && index[2] < _gridSize[2]) {
      newValue = __interpolatePhi(position);
    }
    newPhi[id] = newValue;
  }
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    phi[id] = newPhi[id];
  }
}

void LevelSetFluid3::advectWeno() {}

double LevelSetFluid3::__interpolateVelocityU(Eigen::Array3d position) {
  double _min, _max;
  return __interpolateVelocityU(position, _min, _max);
}

double LevelSetFluid3::__interpolateVelocityU(Eigen::Array3d position,
                                              double &_min, double &_max) {
  Eigen::Array3d h = getH();
  auto &_u = getFaceScalarData(0, _velocityId);
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  // Treating values that are outside domain

  if (position[1] < cellCenter[1])
    index[1]--;
  if (position[2] < cellCenter[2])
    index[2]--;

  // Check tangential values
  if (index[1] < 0)
    index[1]++;
  if (index[1] >= _gridSize[1])
    index[1]--;
  if (index[2] < 0)
    index[2]++;
  if (index[2] >= _gridSize[2])
    index[2]--;

  std::vector<int> iCandidates, jCandidates, kCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
    kCandidates.push_back(index[2] + i);
  }
  std::vector<Eigen::Array3d> points;
  std::vector<double> values;

  for (auto w : kCandidates)
    for (auto v : jCandidates)
      for (auto u : iCandidates) {
        int x, y, z;
        if (v >= 0 && v < _gridSize[1] && w >= 0 && w < _gridSize[2])
          points.emplace_back(u * h[0], (v + 0.5) * h[1], (w + 0.5) * h[2]);
        else {
          int vv = std::max(0, std::min(_gridSize[1] - 1, v));
          int ww = std::max(0, std::min(_gridSize[2] - 1, w));
          if ((v < 0 || v >= _gridSize[1]) && (w < 0 || w >= _gridSize[2]))
            points.emplace_back(u * h[0], vv * h[1], ww * h[2]);
          else if ((v < 0 || v >= _gridSize[1]))
            points.emplace_back(u * h[0], vv * h[1], (w + 0.5) * h[2]);
          else if ((w < 0 || w >= _gridSize[2]))
            points.emplace_back(u * h[0], (v + 0.5) * h[1], ww * h[2]);
        }
        x = std::max(0, std::min(_gridSize[0], u));
        y = std::max(0, std::min(_gridSize[1] - 1, v));
        z = std::max(0, std::min(_gridSize[2] - 1, w));
        values.emplace_back(_u[ijkToid(x, y, z)]);
        _min = std::min(values.back(), _min);
        _max = std::max(values.back(), _max);
      }
  double velocity = Ramuh::Interpolator::tricubic(position, points, values);
  return velocity;
}

double LevelSetFluid3::__interpolateVelocityV(Eigen::Array3d position) {
  double _min, _max;
  return __interpolateVelocityV(position, _min, _max);
}

double LevelSetFluid3::__interpolateVelocityV(Eigen::Array3d position,
                                              double &_min, double &_max) {
  auto &_v = getFaceScalarData(1, _velocityId);
  Eigen::Array3d h = getH();
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  if (position[0] < cellCenter[0])
    index[0]--;
  if (position[2] < cellCenter[2])
    index[2]--;

  // Check tangential values
  if (index[0] < 0)
    index[0]++;
  if (index[0] >= _gridSize[0])
    index[0]--;
  if (index[2] < 0)
    index[2]++;
  if (index[2] >= _gridSize[2])
    index[2]--;

  std::vector<int> iCandidates, jCandidates, kCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
    kCandidates.push_back(index[2] + i);
  }
  std::vector<Eigen::Array3d> points;
  std::vector<double> values;

  for (auto w : kCandidates)
    for (auto v : jCandidates)
      for (auto u : iCandidates) {
        int x, y, z;
        if (u >= 0 && u < _gridSize[0] && w >= 0 && w < _gridSize[2])
          points.emplace_back((u + 0.5) * h[0], v * h[1], (w + 0.5) * h[2]);
        else {
          int uu = std::max(0, std::min(_gridSize[0] - 1, u));
          int ww = std::max(0, std::min(_gridSize[2] - 1, w));
          if ((u < 0 || u >= _gridSize[0]) && (w < 0 || w >= _gridSize[2]))
            points.emplace_back(uu * h[0], v * h[1], ww * h[2]);
          else if ((u < 0 || u >= _gridSize[0]))
            points.emplace_back(uu * h[0], v * h[1], (w + 0.5) * h[2]);
          else if ((w < 0 || w >= _gridSize[2]))
            points.emplace_back((u + 0.5) * h[0], v * h[1], ww * h[2]);
        }
        x = std::max(0, std::min(_gridSize[0] - 1, u));
        y = std::max(0, std::min(_gridSize[1], v));
        z = std::max(0, std::min(_gridSize[2] - 1, w));
        values.emplace_back(_v[ijkToid(x, y, z)]);
        _min = std::min(values.back(), _min);
        _max = std::max(values.back(), _max);
      }
  double velocity = Ramuh::Interpolator::tricubic(position, points, values);
  return velocity;
}

double LevelSetFluid3::__interpolateVelocityW(Eigen::Array3d position) {
  double _min, _max;
  return __interpolateVelocityW(position, _min, _max);
}

double LevelSetFluid3::__interpolateVelocityW(Eigen::Array3d position,
                                              double &_min, double &_max) {
  auto &_w = getFaceScalarData(2, _velocityId);
  Eigen::Array3d h = getH();
  Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
  // Check if inside domain
  Eigen::Array3d cellCenter = index.cast<double>() * h + h / 2.0;
  if (position[0] < cellCenter[0])
    index[0]--;
  if (position[1] < cellCenter[1])
    index[1]--;

  // Check tangential values
  if (index[1] < 0)
    index[1]++;
  if (index[1] >= _gridSize[1])
    index[1]--;
  if (index[0] < 0)
    index[0]++;
  if (index[0] >= _gridSize[0])
    index[0]--;

  std::vector<int> iCandidates, jCandidates, kCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + i);
    jCandidates.push_back(index[1] + i);
    kCandidates.push_back(index[2] + i);
  }
  std::vector<Eigen::Array3d> points;
  std::vector<double> values;

  for (auto w : kCandidates)
    for (auto v : jCandidates)
      for (auto u : iCandidates) {
        int x, y, z;
        if (v >= 0 && v < _gridSize[1] && u >= 0 && u < _gridSize[0])
          points.emplace_back((u + 0.5) * h[0], (v + 0.5) * h[1], w * h[2]);
        else {
          int uu = std::max(0, std::min(_gridSize[0] - 1, u));
          int vv = std::max(0, std::min(_gridSize[1] - 1, v));
          if ((v < 0 || v >= _gridSize[1]) && (u < 0 || u >= _gridSize[0]))
            points.emplace_back(uu * h[0], vv * h[1], w * h[2]);
          else if ((v < 0 || v >= _gridSize[1]))
            points.emplace_back((u + 0.5) * h[0], vv * h[1], w * h[2]);
          else if ((u < 0 || u >= _gridSize[0]))
            points.emplace_back(uu * h[0], (v + 0.5) * h[1], (w + 0.5) * h[2]);
        }
        x = std::max(0, std::min(_gridSize[0] - 1, u));
        y = std::max(0, std::min(_gridSize[1] - 1, v));
        z = std::max(0, std::min(_gridSize[2], w));
        values.emplace_back(_w[ijkToid(x, y, z)]);
        _min = std::min(values.back(), _min);
        _max = std::max(values.back(), _max);
      }
  double velocity = Ramuh::Interpolator::tricubic(position, points, values);
  return velocity;
}

double LevelSetFluid3::__interpolatePhi(Eigen::Array3d position) {
  double _min, _max;
  return __interpolatePhi(position, _min, _max);
}

double LevelSetFluid3::__interpolatePhi(Eigen::Array3d position, double &_min,
                                        double &_max) {
  Eigen::Array3d h = getH();
  auto &_phi = getScalarData(_phiId);
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
        x = std::max(0, std::min(_gridSize[0] - 1, u));
        y = std::max(0, std::min(_gridSize[1] - 1, v));
        z = std::max(0, std::min(_gridSize[2] - 1, w));
        values.emplace_back(_phi[ijkToid(x, y, z)]);

        _min = std::min(values.back(), _min);
        _max = std::max(values.back(), _max);
      }
  double phi = Ramuh::Interpolator::tricubic(position, points, values);
  return std::min(_max, std::max(_min, phi));
  // return phi;
}

} // namespace Leviathan