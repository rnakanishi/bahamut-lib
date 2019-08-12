#include <fluids/levelset_fluid3.h>
#include <blas/interpolator.h>
#include <blas/weno.h>

namespace Leviathan {

LevelSetFluid3::LevelSetFluid3() : Ramuh::MacGrid3() {
  _phiId = newScalarLabel("phi");
  _velocityId = newFaceScalarLabel("velocity");
  //  --> Even though velocities are vector, MAC grid split them into scalar
  //  pieces for each coordinate
  _dt = 1. / 60;
  _tolerance = 1e-10;

  newScalarLabel("newPhi");
}

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
  auto &newPhi = getScalarVector("newPhi");
  auto &phi = getScalarVector(_phiId);

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

void LevelSetFluid3::advectWeno() {
  size_t nCells = cellCount();
  std::vector<double> values(6); // function values
  auto h = getH();

  // For each coordinate, perform HJ-WENO advection separatedly

  auto &_u = getFaceScalarVector(0, _velocityId);
  auto &_v = getFaceScalarVector(1, _velocityId);
  auto &_w = getFaceScalarVector(2, _velocityId);
  auto &_phi = getScalarVector(_phiId);
  auto &newPhi = getScalarVector("newPhi");
#pragma omp parallel for
  for (size_t id = 0; id < nCells; id++) {
    Eigen::Vector3d velocity;
    Eigen::Vector3d dPhi;
    Eigen::Array3i ijk;
    std::tie(ijk[0], ijk[1], ijk[2]) = idToijk(id);

    // Average velocities from faces
    velocity[0] = (_u[id] + _u[ijkToid(ijk[0] + 1, ijk[1], ijk[2])]) / 2;
    velocity[1] = (_v[id] + _v[ijkToid(ijk[0], ijk[1] + 1, ijk[2])]) / 2;
    velocity[2] = (_w[id] + _w[ijkToid(ijk[0], ijk[1], ijk[2] + 1)]) / 2;

    // For center points, check if velocity is positive or negative. If
    // negative, the perform right-weno, otherwise, perform legt-weno.
    // TODO: treat boundary index
    bool isNegative = true;
    if (velocity[0] >= 0)
      isNegative = false;
    values[0] = (isNegative) ? _u[ijkToid(ijk[0] - 2, ijk[1], ijk[2])]
                             : _u[ijkToid(ijk[0] - 3, ijk[1], ijk[2])];
    values[1] = (isNegative) ? _u[ijkToid(ijk[0] - 1, ijk[1], ijk[2])]
                             : _u[ijkToid(ijk[0] - 2, ijk[1], ijk[2])];
    values[2] = (isNegative) ? _u[ijkToid(ijk[0] - 0, ijk[1], ijk[2])]
                             : _u[ijkToid(ijk[0] - 1, ijk[1], ijk[2])];
    values[3] = (isNegative) ? _u[ijkToid(ijk[0] + 1, ijk[1], ijk[2])]
                             : _u[ijkToid(ijk[0] - 0, ijk[1], ijk[2])];
    values[4] = (isNegative) ? _u[ijkToid(ijk[0] + 2, ijk[1], ijk[2])]
                             : _u[ijkToid(ijk[0] + 1, ijk[1], ijk[2])];
    values[5] = (isNegative) ? _u[ijkToid(ijk[0] + 3, ijk[1], ijk[2])]
                             : _u[ijkToid(ijk[0] + 2, ijk[1], ijk[2])];
    dPhi[0] = Ramuh::Weno::evaluate(values, h[0], isNegative);

    isNegative = true;
    if (velocity[1] >= 0)
      isNegative = false;
    values[0] = (isNegative) ? _v[ijkToid(ijk[0], ijk[1] - 2, ijk[2])]
                             : _v[ijkToid(ijk[0], ijk[1] - 3, ijk[2])];
    values[1] = (isNegative) ? _v[ijkToid(ijk[0], ijk[1] - 1, ijk[2])]
                             : _v[ijkToid(ijk[0], ijk[1] - 2, ijk[2])];
    values[2] = (isNegative) ? _v[ijkToid(ijk[0], ijk[1] - 0, ijk[2])]
                             : _v[ijkToid(ijk[0], ijk[1] - 1, ijk[2])];
    values[3] = (isNegative) ? _v[ijkToid(ijk[0], ijk[1] + 1, ijk[2])]
                             : _v[ijkToid(ijk[0], ijk[1] - 0, ijk[2])];
    values[4] = (isNegative) ? _v[ijkToid(ijk[0], ijk[1] + 2, ijk[2])]
                             : _v[ijkToid(ijk[0], ijk[1] + 1, ijk[2])];
    values[5] = (isNegative) ? _v[ijkToid(ijk[0], ijk[1] + 3, ijk[2])]
                             : _v[ijkToid(ijk[0], ijk[1] + 2, ijk[2])];
    dPhi[1] = Ramuh::Weno::evaluate(values, h[0], isNegative);

    isNegative = true;
    if (velocity[2] >= 0)
      isNegative = false;
    values[0] = (isNegative) ? _w[ijkToid(ijk[0], ijk[1], ijk[2] - 2)]
                             : _w[ijkToid(ijk[0], ijk[1], ijk[2] - 3)];
    values[1] = (isNegative) ? _w[ijkToid(ijk[0], ijk[1], ijk[2] - 1)]
                             : _w[ijkToid(ijk[0], ijk[1], ijk[2] - 2)];
    values[2] = (isNegative) ? _w[ijkToid(ijk[0], ijk[1], ijk[2] - 0)]
                             : _w[ijkToid(ijk[0], ijk[1], ijk[2] - 1)];
    values[3] = (isNegative) ? _w[ijkToid(ijk[0], ijk[1], ijk[2] + 1)]
                             : _w[ijkToid(ijk[0], ijk[1], ijk[2] - 0)];
    values[4] = (isNegative) ? _w[ijkToid(ijk[0], ijk[1], ijk[2] + 2)]
                             : _w[ijkToid(ijk[0], ijk[1], ijk[2] + 1)];
    values[5] = (isNegative) ? _w[ijkToid(ijk[0], ijk[1], ijk[2] + 3)]
                             : _w[ijkToid(ijk[0], ijk[1], ijk[2] + 2)];
    dPhi[2] = Ramuh::Weno::evaluate(values, h[0], isNegative);

    // Proceed to time integration and material derivative
    // Euler method
    newPhi[id] = -velocity.dot(dPhi) * _dt;
    newPhi[id] = (newPhi[id] + _phi[id]);
  }

#pragma omp parallel for
  for (size_t id = 0; id < nCells; id++) {
    _phi[id] = newPhi[id];
  }
}

double LevelSetFluid3::__interpolateVelocityU(Eigen::Array3d position) {
  double _min, _max;
  __interpolateVelocityU(position, _min, _max);
}

double LevelSetFluid3::__interpolateVelocityU(Eigen::Array3d position,
                                              double &_min, double &_max) {
  Eigen::Array3d h = getH();
  auto &_u = getFaceScalarVector(0, _velocityId);
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
  __interpolateVelocityV(position, _min, _max);
}

double LevelSetFluid3::__interpolateVelocityV(Eigen::Array3d position,
                                              double &_min, double &_max) {
  auto &_v = getFaceScalarVector(1, _velocityId);
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
  __interpolateVelocityW(position, _min, _max);
}

double LevelSetFluid3::__interpolateVelocityW(Eigen::Array3d position,
                                              double &_min, double &_max) {
  auto &_w = getFaceScalarVector(2, _velocityId);
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
  auto &_phi = getScalarVector(_phiId);
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