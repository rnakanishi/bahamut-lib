#include <fluids/levelset_fluid3.h>
#include <blas/interpolator.h>
#include <blas/weno.h>
#include <cmath>

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
  // FIXME: Possibly contain errors. Track and fix them
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

void LevelSetFluid3::advectUpwind() {
  auto h = getH();
  auto &phi = getScalarData(_phiId);
  std::vector<double> upwind(cellCount());

#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    std::vector<double> values(6);
    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    double dotGrad = 0.;

    for (size_t coord = 0; coord < 3; coord++) {
      auto &u = getFaceScalarData(coord, _velocityId);
      auto faceijk = faceIdToijk(coord, id);
      int facei = faceijk[0], facej = faceijk[1], facek = faceijk[2];
      double velocity;

      int cellijk;
      if (coord == 0)
        cellijk = i;
      else if (coord == 1)
        cellijk = j;
      else if (coord == 2)
        cellijk = k;

      // For each coordinate, compute its own velocity
      switch (coord) {
      case 0: // X coordinate
        velocity = ((u[faceijkToid(coord, facei, facej, facek)] +
                     u[faceijkToid(coord, facei + 1, facej, facek)]) /
                    2);
        break;
      case 1: // Y coordinate
        velocity = ((u[faceijkToid(coord, facei, facej, facek)] +
                     u[faceijkToid(coord, facei, facej + 1, facek)]) /
                    2);
        break;
      case 2: // Z coordinate
        velocity = ((u[faceijkToid(coord, facei, facej, facek)] +
                     u[faceijkToid(coord, facei, facej, facek + 1)]) /
                    2);
        break;
      }
      if (velocity <= 0) {
        int index;
        switch (coord) {
        case 0:
          index = std::max(0, i - 1);
          dotGrad +=
              velocity * (phi[id] - phi[ijkToid(index, j, k)]) / h[coord];
          break;
        case 1:
          index = std::max(0, j - 1);
          dotGrad +=
              velocity * (phi[id] - phi[ijkToid(i, index, k)]) / h[coord];
          break;
        case 2:
          index = std::max(0, k - 1);
          dotGrad +=
              velocity * (phi[id] - phi[ijkToid(i, j, index)]) / h[coord];
          break;
        }
      } else {
        int index;
        switch (coord) {
        case 0:
          index = std::min(_gridSize[0] - 1, i + 1);
          dotGrad +=
              velocity * (phi[ijkToid(index, j, k)] - phi[id]) / h[coord];
          break;
        case 1:
          index = std::min(_gridSize[1] - 1, j + 1);
          dotGrad +=
              velocity * (phi[ijkToid(i, index, k)] - phi[id]) / h[coord];
          break;
        case 2:
          index = std::min(_gridSize[2] - 1, k + 1);
          dotGrad +=
              velocity * (phi[ijkToid(i, j, index)] - phi[id]) / h[coord];
          break;
        }
      }
    }
    upwind[id] = dotGrad;
  }

// Time integration step. Euler method
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    phi[id] = phi[id] - upwind[id] * _dt;
  }
}

void LevelSetFluid3::advectWeno() {
  auto h = getH();
  auto &phi = getScalarData(_phiId);
  std::vector<double> weno(cellCount());

// Weno computation
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    std::vector<double> values(6);
    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    weno[id] = 0;

    for (size_t coord = 0; coord < 3; coord++) {
      auto &u = getFaceScalarData(coord, _velocityId);
      auto faceijk = ijk; // FIXME: change to faceIdToijk(id) and check errors
      int facei = faceijk[0], facej = faceijk[1], facek = faceijk[2];
      double velocity;

      int cellijk;
      if (coord == 0)
        cellijk = i;
      else if (coord == 1)
        cellijk = j;
      else if (coord == 2)
        cellijk = k;

      // For each coordinate, compute its own velocity
      switch (coord) {
      case 0: // X coordinate
        velocity = ((u[faceijkToid(coord, facei, facej, facek)] +
                     u[faceijkToid(coord, facei + 1, facej, facek)]) /
                    2);
        break;
      case 1: // Y coordinate
        velocity = ((u[faceijkToid(coord, facei, facej, facek)] +
                     u[faceijkToid(coord, facei, facej + 1, facek)]) /
                    2);
        break;
      case 2: // Z coordinate
        velocity = ((u[faceijkToid(coord, facei, facej, facek)] +
                     u[faceijkToid(coord, facei, facej, facek + 1)]) /
                    2);
        break;
      }
      // For positive velocity  upwind-based scheme is used. Otherwise,
      // downwind
      if (velocity <= 0) {
        // ========= DOWNWIND =========
        for (int ival = 0, inc = 3; ival < 7; ival++, inc--) {
          int index, index1, faceIndex;
          // Each coordinate have to be treated separatedly due to cell-face
          // indexation
          switch (coord) {
            // FIXME: for each cell in the stencil, use average velocity,
            // instead of face velocity
          case 0:
            index = std::min(_gridSize[coord] - 1, std::max(1, i + inc));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, i + inc - 1));
            faceIndex = std::min(_gridSize[coord], std::max(0, facei + inc));
            values[ival] =
                u[faceijkToid(coord, faceIndex, facej, facek)] *
                (phi[ijkToid(index, j, k)] - phi[ijkToid(index1, j, k)]) /
                h[coord];
            break;
          case 1:
            index = std::min(_gridSize[coord] - 1, std::max(1, j + inc));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, j + inc - 1));
            faceIndex = std::min(_gridSize[coord], std::max(0, facej + inc));
            values[ival] =
                u[faceijkToid(coord, facei, faceIndex, facek)] *
                (phi[ijkToid(i, index, k)] - phi[ijkToid(i, index1, k)]) /
                h[coord];
            break;
          case 2:
            index = std::min(_gridSize[coord] - 1, std::max(1, k + inc));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, k + inc - 1));
            faceIndex = std::min(_gridSize[coord], std::max(0, facek + inc));
            values[ival] =
                u[faceijkToid(coord, facei, facej, faceIndex)] *
                (phi[ijkToid(i, j, index)] - phi[ijkToid(i, j, index1)]) /
                h[coord];
            break;
          }
        }
        // Boundary conditions. This ensure weno doesnt use outside values
        if (cellijk <= 1)
          values[5] = 1e4;
        if (cellijk <= 0)
          values[4] = 1e4;
        if (cellijk >= _gridSize[coord] - 2)
          values[0] = 1e4;
        if (cellijk >= _gridSize[coord] - 1)
          values[1] = 1e4;
      } else {
        // ========= UPWIND =========
        for (int ival = 0, inc = -3; ival < 7; ival++, inc++) {
          int index, index1, faceIndex;
          // Each coordinate have to be treated separatedly due to cell-face
          // indexation
          switch (coord) {
          case 0:
            index = std::min(_gridSize[coord] - 1, std::max(1, i + inc + 1));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, i + inc));
            faceIndex =
                std::min(_gridSize[coord], std::max(0, facei + inc + 1));
            values[ival] =
                u[faceijkToid(coord, faceIndex, facej, facek)] *
                (phi[ijkToid(index, j, k)] - phi[ijkToid(index1, j, k)]) /
                h[coord];
            break;
          case 1:
            index = std::min(_gridSize[coord] - 1, std::max(1, j + inc + 1));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, j + inc));
            faceIndex =
                std::min(_gridSize[coord], std::max(0, facej + inc + 1));
            values[ival] =
                u[faceijkToid(coord, facei, faceIndex, facek)] *
                (phi[ijkToid(i, index, k)] - phi[ijkToid(i, index1, k)]) /
                h[coord];
            break;
          case 2:
            index = std::min(_gridSize[coord] - 1, std::max(1, k + inc + 1));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, k + inc));
            faceIndex =
                std::min(_gridSize[coord], std::max(0, facek + inc + 1));
            values[ival] =
                u[faceijkToid(coord, facei, facej, faceIndex)] *
                (phi[ijkToid(i, j, index)] - phi[ijkToid(i, j, index1)]) /
                h[coord];
            break;
          }
        }
        if (cellijk <= 1)
          values[0] = 1e4;
        if (cellijk <= 0)
          values[1] = 1e4;
        if (cellijk >= _gridSize[coord] - 2)
          values[5] = 1e4;
        if (cellijk >= _gridSize[coord] - 1)
          values[4] = 1e4;
      }

      double dPhi = Ramuh::Weno::evaluate(values, h[coord], false);
      // Summing up all gradients to obtain flux
      weno[id] += dPhi;
    }
  }

// Time integration step. Euler method
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    phi[id] = phi[id] - weno[id] * _dt;
  }
}

void LevelSetFluid3::redistance() {
  auto &phi = getScalarData(_phiId);
  std::vector<double> gradient(cellCount());
  std::vector<double> smooth(cellCount());
  std::vector<double> phiSign(cellCount());
  double eps = 1e-16;
  for (size_t id = 0; id < cellCount(); id++)
    phiSign[id] = (phi[id] >= 0) ? 1 : -1;

  for (size_t t = 0; t < 20; t++) {
    // FIXME: Run this equaiton until stationary state is reached
    // #pragma omp parallel for
    for (size_t id = 0; id < cellCount(); id++) {
      // Smooth function
      smooth[id] = phi[id] / sqrt(phi[id] * phi[id] + eps * eps);

      // Compute gradients
      std::vector<double> partials(3);
      auto h = getH();
      auto ijk = idToijk(id);
      int i = ijk[0], j = ijk[1], k = ijk[2];
      gradient[id] = 0.0;
      int index;
      if (phiSign[id] > 0) {
        partials[0] = (phi[id] - phi[ijkToid(std::max(0, i - 1), j, k)]) / h[0];
        partials[1] = (phi[id] - phi[ijkToid(i, std::max(0, j - 1), k)]) / h[1];
        partials[2] = (phi[id] - phi[ijkToid(i, j, std::max(0, k - 1))]) / h[2];
      } else {
        partials[0] =
            (phi[id] - phi[ijkToid(std::min(_gridSize[0] - 1, i + 1), j, k)]) /
            h[0];
        partials[1] =
            (phi[ijkToid(i, std::min(_gridSize[1] - 1, j + 1), k)] - phi[id]) /
            h[1];
        partials[2] =
            (phi[ijkToid(i, j, std::min(_gridSize[2] - 1, k + 1))] - phi[id]) /
            h[2];
      }

      if (abs(phi[id]) >= 1e-10) {
        gradient[id] =
            sqrt(partials[0] * partials[0] + partials[1] * partials[1] +
                 partials[2] * partials[2]) -
            1;
      }
    }

    for (size_t id = 0; id < cellCount(); id++) {
      phi[id] = phi[id] - _dt * smooth[id] * gradient[id];
    }
  }
}

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