#include <fluids/levelset_fluid3.h>
#include <blas/interpolator.h>
#include <blas/weno.h>
#include <cmath>
#include <iostream>

namespace Leviathan {

LevelSetFluid3::LevelSetFluid3()
    : LevelSetFluid3(Eigen::Array3i(32, 32, 32),
                     Ramuh::BoundingBox3(Eigen::Array3d(-1, -1, -1),
                                         Eigen::Array3d(1, 1, 1))) {}

LevelSetFluid3::LevelSetFluid3(Eigen::Array3i gridSize,
                               Ramuh::BoundingBox3 domain)
    : Ramuh::MacGrid3(domain, gridSize) {
  _phiId = newCellScalarLabel("phi", 1e8);
  _cellVelocityId = newCellScalarLabel("cellVelocity");
  _faceVelocityId = newFaceScalarLabel("faceVelocity");
  //  --> Even though velocities are vector, MAC grid split them into scalar
  //  pieces for each coordinate
  _originalDt = _dt = 1. / 60;
  _ellapsedDt = 0;
  _tolerance = 1e-10;

  _cellGradientId = newCellArrayLabel("cellGradient");
  newCellScalarLabel("newPhi");
}

void LevelSetFluid3::computeCellsGradient() {
  auto &phi = getCellScalarData("phi");
  auto &gradient = getCellArrayData("cellGradient");
  auto h = getH();
#pragma omp parallel for
  for (int k = 0; k < _gridSize[2]; k++)
    for (int j = 0; j < _gridSize[1]; j++)
      for (int i = 0; i < _gridSize[0]; i++) {
        int id = ijkToid(i, j, k);
        if (i < _gridSize[0])
          gradient[id][0] = (phi[ijkToid(i + 1, j, k)] - phi[id]) / h[0];
        else
          gradient[id] = (phi[id] - phi[ijkToid(i - 1, j, k)]) / h[0];

        if (j < _gridSize[1])
          gradient[id][1] = (phi[ijkToid(i, j + 1, k)] - phi[id]) / h[1];
        else
          gradient[id] = (phi[id] - phi[ijkToid(i, j - 1, k)]) / h[1];

        if (k < _gridSize[2])
          gradient[id][2] = (phi[ijkToid(i, j, k + 1)] - phi[id]) / h[2];
        else
          gradient[id] = (phi[id] - phi[ijkToid(i, j, k - 1)]) / h[2];
      }
}

void LevelSetFluid3::advectSemiLagrangean() {
  size_t nCells = cellCount();
  auto &newPhi = getCellScalarData("newPhi");
  auto &phi = getCellScalarData(_phiId);

#pragma omp parallel for
  for (size_t id = 0; id < nCells; id++) {
    Eigen::Array3d position, h, velocity;
    h = getH();

    velocity[0] = interpolateFaceScalarData(0, _faceVelocityId, position);
    velocity[1] = interpolateFaceScalarData(1, _faceVelocityId, position);
    velocity[2] = interpolateFaceScalarData(2, _faceVelocityId, position);

    position = position - (velocity * _dt);
    Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
    // Check if inside domain
    double newValue = phi[id];
    if ((velocity.matrix().norm() > _tolerance) && index[0] >= 0 &&
        index[1] >= 0 && index[2] >= 0 && index[0] < _gridSize[0] &&
        index[1] < _gridSize[1] && index[2] < _gridSize[2]) {
      newValue = interpolateCellScalarData(_phiId, position);
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
  auto &phi = getCellScalarData(_phiId);
  std::vector<double> upwind(cellCount());

#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    std::vector<double> values(6);
    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    double dotGrad = 0.;

    for (size_t coord = 0; coord < 3; coord++) {
      auto &uvw = getFaceScalarData(coord, _faceVelocityId);
      auto faceijk = ijk;
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
        velocity = ((uvw[faceijkToid(coord, facei, facej, facek)] +
                     uvw[faceijkToid(coord, facei + 1, facej, facek)]) /
                    2);
        break;
      case 1: // Y coordinate
        velocity = ((uvw[faceijkToid(coord, facei, facej, facek)] +
                     uvw[faceijkToid(coord, facei, facej + 1, facek)]) /
                    2);
        break;
      case 2: // Z coordinate
        velocity = ((uvw[faceijkToid(coord, facei, facej, facek)] +
                     uvw[faceijkToid(coord, facei, facej, facek + 1)]) /
                    2);
        break;
      }
      if (velocity > 0) { // DOWNWIND
        int index;
        switch (coord) {
        case 0:
          if (i > 0)
            dotGrad +=
                velocity * (phi[id] - phi[ijkToid(i - 1, j, k)]) / h[coord];
          else
            dotGrad +=
                velocity * (phi[ijkToid(i + 1, j, k) - phi[id]]) / h[coord];
          break;
        case 1:
          if (j > 0)
            dotGrad +=
                velocity * (phi[id] - phi[ijkToid(i, j - 1, k)]) / h[coord];
          else
            dotGrad +=
                velocity * (phi[ijkToid(i, j + 1, k)] - phi[id]) / h[coord];
          break;
        case 2:
          if (k > 0)
            dotGrad +=
                velocity * (phi[id] - phi[ijkToid(i, j, k - 1)]) / h[coord];
          else
            dotGrad +=
                velocity * (phi[ijkToid(i, j, k + 1)] - phi[id]) / h[coord];
          break;
        }
      } else { // UP
        int index;
        switch (coord) {
        case 0:
          if (i < _gridSize[0] - 1)
            dotGrad +=
                velocity * (phi[ijkToid(i + 1, j, k)] - phi[id]) / h[coord];
          else
            dotGrad +=
                velocity * (phi[id] - phi[ijkToid(i - 1, j, k)]) / h[coord];
          break;
        case 1:
          if (j < _gridSize[1] - 1)
            dotGrad +=
                velocity * (phi[ijkToid(i, j + 1, k)] - phi[id]) / h[coord];
          else
            dotGrad +=
                velocity * (phi[id] - phi[ijkToid(i, j - 1, k)]) / h[coord];
          break;
        case 2:
          if (k < _gridSize[2] - 1)
            dotGrad +=
                velocity * (phi[ijkToid(i, j, k + 1)] - phi[id]) / h[coord];
          else
            dotGrad +=
                velocity * (phi[id] - phi[ijkToid(i, j, k - 1)]) / h[coord];
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

void LevelSetFluid3::computeWenoGradient() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  auto &gradient = getCellScalarData(_cellGradientId);

// Weno computation
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    std::vector<double> values(6);
    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    gradient[id] = 0;

    for (size_t coord = 0; coord < 3; coord++) {
      auto &uvw = getFaceScalarData(coord, _faceVelocityId);
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
        velocity = ((uvw[faceijkToid(coord, facei, facej, facek)] +
                     uvw[faceijkToid(coord, facei + 1, facej, facek)]) /
                    2);
        break;
      case 1: // Y coordinate
        velocity = ((uvw[faceijkToid(coord, facei, facej, facek)] +
                     uvw[faceijkToid(coord, facei, facej + 1, facek)]) /
                    2);
        break;
      case 2: // Z coordinate
        velocity = ((uvw[faceijkToid(coord, facei, facej, facek)] +
                     uvw[faceijkToid(coord, facei, facej, facek + 1)]) /
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
                (phi[ijkToid(index, j, k)] - phi[ijkToid(index1, j, k)]) /
                h[coord];
            break;
          case 1:
            index = std::min(_gridSize[coord] - 1, std::max(1, j + inc));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, j + inc - 1));
            faceIndex = std::min(_gridSize[coord], std::max(0, facej + inc));
            values[ival] =
                (phi[ijkToid(i, index, k)] - phi[ijkToid(i, index1, k)]) /
                h[coord];
            break;
          case 2:
            index = std::min(_gridSize[coord] - 1, std::max(1, k + inc));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, k + inc - 1));
            faceIndex = std::min(_gridSize[coord], std::max(0, facek + inc));
            values[ival] =
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
                (phi[ijkToid(index, j, k)] - phi[ijkToid(index1, j, k)]) /
                h[coord];
            break;
          case 1:
            index = std::min(_gridSize[coord] - 1, std::max(1, j + inc + 1));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, j + inc));
            faceIndex =
                std::min(_gridSize[coord], std::max(0, facej + inc + 1));
            values[ival] =
                (phi[ijkToid(i, index, k)] - phi[ijkToid(i, index1, k)]) /
                h[coord];
            break;
          case 2:
            index = std::min(_gridSize[coord] - 1, std::max(1, k + inc + 1));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, k + inc));
            faceIndex =
                std::min(_gridSize[coord], std::max(0, facek + inc + 1));
            values[ival] =
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
      gradient[id] += dPhi;
    }
  }
}

void LevelSetFluid3::advectWeno() {
  auto &gradient = getCellArrayData(_cellGradientId);
  auto &cellVelocity = getCellArrayData(_cellVelocityId);
  auto &phi = getCellScalarData(_phiId);

  std::vector<double> phiN(phi.size());
  for (size_t i = 0; i < phi.size(); i++) {
    phiN[i] = phi[i];
  }

  // computeCentralGradient();
  computeCellVelocity();
  computeWenoGradient();
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    phi[id] =
        phi[id] - cellVelocity[id].matrix().dot(gradient[id].matrix()) * _dt;
  }

  // computeCentralGradient();
  computeCellVelocity();
  computeWenoGradient();
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    phi[id] =
        phi[id] - cellVelocity[id].matrix().dot(gradient[id].matrix()) * _dt;
  }

  for (size_t i = 0; i < phi.size(); i++) {
    phi[i] = 0.75 * phiN[i] + 0.25 * phi[i];
  }

  // computeCentralGradient();
  computeCellVelocity();
  computeWenoGradient();
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    phi[id] =
        phi[id] - cellVelocity[id].matrix().dot(gradient[id].matrix()) * _dt;
  }

#pragma omp parallel for
  for (size_t i = 0; i < phi.size(); i++) {
    phi[i] = phiN[i] / 3 + 2 * phi[i] / 3;
  }
}

void LevelSetFluid3::computeCellVelocity() {
  auto &cellVelocity = getCellArrayData(_cellVelocityId);

#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];

    for (size_t coord = 0; coord < 2; coord++) {
      auto &uvw = getFaceScalarData(coord, _faceVelocityId);
      double velocity;
      switch (coord) {
      case 0:
        velocity = (uvw[faceijkToid(coord, i, j, k)] +
                    uvw[faceijkToid(coord, i + 1, j, k)]) /
                   2;
        break;
      case 1:
        velocity = (uvw[faceijkToid(coord, i, j, k)] +
                    uvw[faceijkToid(coord, i, j + 1, k)]) /
                   2;
        break;
      case 2:
        velocity = (uvw[faceijkToid(coord, i, j, k)] +
                    uvw[faceijkToid(coord, i, j, k + 1)]) /
                   2;
        break;
      }
      cellVelocity[id][coord] = velocity;
    }
  }
}

void LevelSetFluid3::advectCip() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  auto &G = getCellArrayData(_cellGradientId);
  auto &uVelocity = getFaceScalarData(0, _cellVelocityId);
  auto &vVelocity = getFaceScalarData(1, _cellVelocityId);
  auto &wVelocity = getFaceScalarData(2, _cellVelocityId);
  static int count = 0;
  std::vector<double> newPhi(cellCount(), 0);
  std::vector<Eigen::Array3d> newGrad(cellCount(), Eigen::Array3d(0));

  // #pragma omp parallel for
  for (size_t id = 0; id < cellCount(); id++) {
    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];

    double uCellVelocity = 0.5 * (uVelocity[faceijkToid(0, i, j, k)] +
                                  uVelocity[faceijkToid(0, i + 1, j, k)]);
    double vCellVelocity = 0.5 * (vVelocity[faceijkToid(1, i, j, k)] +
                                  vVelocity[faceijkToid(1, i, j + 1, k)]);
    double wCellVelocity = 0.5 * (wVelocity[faceijkToid(2, i, j, k)] +
                                  wVelocity[faceijkToid(2, i, j, k + 1)]);

    double XX = -uCellVelocity * _dt;
    double YY = -vCellVelocity * _dt;
    double ZZ = -wCellVelocity * _dt;
    int isgn = (uCellVelocity >= 0) ? 1 : -1;
    int jsgn = (vCellVelocity >= 0) ? 1 : -1;
    int ksgn = (wCellVelocity >= 0) ? 1 : -1;
    double dx = h[0] * isgn;
    double dy = h[1] * jsgn;
    double dz = h[2] * ksgn;
    int im1 = std::max(0, std::min(_gridSize[0] - 1, i - isgn));
    int jm1 = std::max(0, std::min(_gridSize[1] - 1, j - jsgn));
    int km1 = std::max(0, std::min(_gridSize[2] - 1, k - ksgn));
    int idim1 = ijkToid(im1, j, k);
    int idjm1 = ijkToid(i, jm1, k);
    int idkm1 = ijkToid(i, j, km1);

    double B[19];
    B[16] = -phi[id] + phi[idim1] + phi[idjm1] - phi[ijkToid(im1, jm1, k)];
    B[17] = -phi[id] + phi[idim1] + phi[idkm1] - phi[ijkToid(im1, j, km1)];
    B[18] = -phi[id] + phi[idjm1] + phi[idkm1] - phi[ijkToid(i, jm1, km1)];
    B[0] = (-2 * (phi[idim1] - phi[id]) + (G[idim1][0] + G[id][0]) * dx) /
           (pow(dx, 3));
    B[1] = -(B[16] + (G[idjm1][0] - G[id][0]) * dx) / (dx * dx * dy);
    B[2] = -(B[17] + (G[idkm1][0] - G[id][0]) * dx) / (dx * dx * dz);
    B[3] = (3 * (phi[idim1] - phi[id]) - (G[idim1][0] + 2 * G[id][0]) * dx) /
           (dx * dx);
    B[4] = (B[16] + (G[idjm1][0] - G[id][0]) * dx +
            (G[idim1][1] - G[id][1]) * dy) /
           (dx * dy);
    B[5] = (-2 * (phi[idjm1] - phi[id]) + (G[idjm1][1] + G[id][1]) * dy) /
           (dy * dy * dy);
    B[6] = -(B[18] + (G[idkm1][1] - G[id][1]) * dy) / (dy * dy * dz);
    B[7] = -(B[16] + (G[idim1][1] - G[id][1]) * dy) / (dx * dy * dy);
    B[8] = (3 * (phi[idjm1] - phi[id]) - (G[idim1][1] + 2 * G[id][1]) * dy) /
           (dy * dy);
    B[9] = (B[18] + (G[idkm1][1] - G[id][1]) * dy +
            (G[idjm1][2] - G[id][2]) * dz) /
           (dy * dz);
    B[10] = (-2 * (phi[idkm1] - phi[id]) + (G[idkm1][2] + G[id][2]) * dz) /
            (dz * dz * dz);
    B[11] = -(B[17] + (G[idim1][2] - G[id][2]) * dz) / (dx * dz * dz);
    B[12] = -(B[18] + (G[idjm1][2] - G[id][2]) * dz) / (dy * dz * dz);
    B[13] = (3 * (phi[idkm1] - phi[id]) - (G[idkm1][2] + 2 * G[id][2]) * dz) /
            (dz * dz);
    B[14] = (B[17] + (G[idim1][2] - G[id][2]) * dz +
             (G[idkm1][0] - G[id][0]) * dx) /
            (dx * dz);
    B[15] = (B[16] + phi[idkm1] - phi[ijkToid(i, jm1, km1)] -
             phi[ijkToid(im1, j, km1)] + phi[ijkToid(im1, jm1, km1)]) /
            (dx * dy * dz);

    newPhi[id] = 0;
    newPhi[id] += XX * ((B[0] * XX + B[1] * YY + B[2] * ZZ + B[3]) * XX +
                        B[4] * YY + G[id][0]);
    newPhi[id] += YY * ((B[5] * YY + B[6] * ZZ + B[7] * XX + B[8]) * YY +
                        B[9] * ZZ + G[id][1]);
    newPhi[id] += ZZ * ((B[10] * ZZ + B[11] * XX + B[12] * YY + B[13]) * ZZ +
                        B[14] * XX + G[id][2]);
    newPhi[id] += B[15] * XX * YY * ZZ + phi[id];
    newGrad[id][0] =
        XX * (3 * XX * B[0] + 2 * YY * B[1] + 2 * ZZ * B[2] + 2 * B[3]) +
        YY * (B[4] + YY * B[7]) + ZZ * (ZZ * B[11] + B[14]) + YY * ZZ * B[15] +
        G[id][0];
    newGrad[id][1] =
        XX * (3 * YY * B[5] + 2 * ZZ * B[6] + 2 * XX * B[7] + 2 * B[8]) +
        XX * (B[4] + XX * B[1]) + ZZ * (ZZ * B[12] + B[9]) + XX * ZZ * B[15] +
        G[id][1];
    newGrad[id][2] =
        ZZ * (3 * ZZ * B[10] + 2 * XX * B[11] + 2 * YY * B[12] + 2 * B[13]) +
        XX * (B[14] + XX * B[2]) + YY * (YY * B[6] + B[9]) + XX * ZZ * B[15] +
        G[id][2];
  }
  for (size_t id = 0; id < cellCount(); id++) {
    phi[id] = newPhi[id];
    G[id] = newGrad[id];
  }
}

void LevelSetFluid3::redistance() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  double eps = h[0];
  double dt = 0.5 * h[0];
  std::vector<double> initialPhi(cellCount());
  std::vector<double> cellSignal(cellCount());
  std::vector<double> gradient(cellCount());
  std::vector<double> interfaceFactor(cellCount());
  std::vector<bool> isInterface(cellCount(), false);
  std::vector<Eigen::Vector3d> direction(cellCount());

#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    initialPhi[id] = phi[id];
    if (phi[id] > 0)
      cellSignal[id] = 1;
    else if (phi[id] < 0)
      cellSignal[id] = -1;
    else
      cellSignal[id] = 0;

    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    interfaceFactor[id] = 0;
    if ((i > 0 && i < _gridSize[0] - 1) && (j > 0 && j < _gridSize[1] - 1) &&
        (k > 0 && k < _gridSize[2] - 1)) {
      if (phi[id] * phi[ijkToid(i - 1, j, k)] <= 0 ||
          phi[id] * phi[ijkToid(i + 1, j, k)] <= 0 ||
          phi[id] * phi[ijkToid(i, j - 1, k)] <= 0 ||
          phi[id] * phi[ijkToid(i, j + 1, k)] <= 0 ||
          phi[id] * phi[ijkToid(i, j, k - 1)] <= 0 ||
          phi[id] * phi[ijkToid(i, j, k + 1)] <= 0) {
        isInterface[id] = true;
        interfaceFactor[id] = 2 * h[0] * phi[id];
        double dx, dy, dz;
        dx = phi[ijkToid(i + 1, j, k)] - phi[ijkToid(i - 1, j, k)];
        dy = phi[ijkToid(i, j + 1, k)] - phi[ijkToid(i, j - 1, k)];
        dz = phi[ijkToid(i, j, k + 1)] - phi[ijkToid(i, j, k - 1)];
        interfaceFactor[id] /= sqrt(dx * dx + dy * dy + dz * dz);
      }
    }
  }
  int t;
  double error = 1;
  for (t = 0; t < 500; t++) {
#pragma omp parallel for
    for (int id = 0; id < cellCount(); id++) {
      auto ijk = idToijk(id);
      int i = ijk[0], j = ijk[1], k = ijk[2];
      double dx[2], dy[2], dz[2]; // Index 0: Dx-, Index 1: Dx+
      dx[0] = dx[1] = dy[0] = dy[1] = dz[0] = dz[1] = 0.;
      {
        if (i > 0)
          dx[0] = (phi[id] - phi[ijkToid(i - 1, j, k)]) / h[0];
        if (i < _gridSize[0] - 1)
          dx[1] = (phi[ijkToid(i + 1, j, k)] - phi[id]) / h[0];

        if (j > 0)
          dy[0] = (phi[id] - phi[ijkToid(i, j - 1, k)]) / h[1];
        if (j < _gridSize[1] - 1)
          dy[1] = (phi[ijkToid(i, j + 1, k)] - phi[id]) / h[1];

        if (k > 0)
          dz[0] = (phi[id] - phi[ijkToid(i, j, k - 1)]) / h[2];
        if (k < _gridSize[2] - 1)
          dz[1] = (phi[ijkToid(i, j, k + 1)] - phi[id]) / h[2];
      }

      double gx, gy, gz;
      if (initialPhi[id] > 0) {
        double a, b;
        a = std::max(0., dx[0]);
        b = std::min(0., dx[1]);
        gx = std::max(a * a, b * b);
        a = std::max(0., dy[0]);
        b = std::min(0., dy[1]);
        gy = std::max(a * a, b * b);
        a = std::max(0., dz[0]);
        b = std::min(0., dz[1]);
        gz = std::max(a * a, b * b);
      } else if (initialPhi[id] < 0) {
        double a, b;
        a = std::min(0., dx[0]);
        b = std::max(0., dx[1]);
        gx = std::max(a * a, b * b);
        a = std::min(0., dy[0]);
        b = std::max(0., dy[1]);
        gy = std::max(a * a, b * b);
        a = std::min(0., dz[0]);
        b = std::max(0., dz[1]);
        gz = std::max(a * a, b * b);
      } else {
        gx = gy = gz = 1.0 / 3;
      }
      gradient[id] = std::sqrt(gx + gy + gz) - 1;
    }

    error = 0.;
// Time integration step. Euler method
#pragma omp parallel for reduction(+ : error)
    for (int id = 0; id < cellCount(); id++) {
      double newPhi = 0.;

      if (!isInterface[id]) {
        newPhi = phi[id] - dt * cellSignal[id] * gradient[id];
      } else {
        newPhi = phi[id] - (dt / h[0]) * (cellSignal[id] * abs(phi[id]) -
                                          interfaceFactor[id]);
      }
      error += abs(phi[id] - newPhi);
      phi[id] = newPhi;
    }
    if (error / cellCount() < dt * h[0] * h[0])
      break;
  }
  std::cerr << "Run redistance for " << t + 1 << " iterations\n";
}

bool LevelSetFluid3::advanceTime() {
  _ellapsedDt += _dt;
  if (_ellapsedDt < _originalDt)
    return false;
  _ellapsedDt = 0.0;
  _dt = _originalDt;
  return true;
}

void LevelSetFluid3::applyCfl() {
  Eigen::Vector2d maxVel = Eigen::Vector2d(0.0, 0.0);
  auto &u = getFaceScalarData(0, _cellVelocityId);
  auto &v = getFaceScalarData(1, _cellVelocityId);
  auto &w = getFaceScalarData(2, _cellVelocityId);
  auto h = getH();

  for (int id = 0; id < cellCount(); id++) {
    auto ijk = idToijk(id);
    int i, j, k;
    i = ijk[0];
    j = ijk[1];
    k = ijk[2];

    Eigen::Vector2d vel;
    vel[0] = (u[id] + u[faceijkToid(0, i + 1, j, k)]) / 2.0;
    vel[1] = (v[id] + v[faceijkToid(1, i, j + 1, k)]) / 2.0;
    vel[1] = (w[id] + w[faceijkToid(2, i, j, k + 1)]) / 2.0;

    if (maxVel.norm() < vel.norm())
      maxVel = vel;
  }
  // Check if cfl condition applies
  // Half timestep if so
  if (maxVel.norm() * (_originalDt - _ellapsedDt) > 0.8 * h[0]) {
    int pieces = std::floor(maxVel.norm() * _dt / (0.9 * h[0])) + 1;
    _dt = _dt / pieces;
  }
}

std::vector<int> LevelSetFluid3::trackSurface() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  double eps = h[0], error = 1e8, dt = 0.5 * h[0];
  std::vector<double> cellSignal(cellCount());
  std::vector<bool> isInterface(cellCount(), false);
  std::vector<int> surface;

  std::fill(_surfaceCells.begin(), _surfaceCells.end(), false);
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    if (i > 0 && i < _gridSize[0] - 1 && j > 0 && j < _gridSize[1] - 1 &&
        k > 0 && k < _gridSize[2] - 1) {
      if (phi[id] * phi[ijkToid(i - 1, j, k)] <= 0 ||
          phi[id] * phi[ijkToid(i + 1, j, k)] <= 0 ||
          phi[id] * phi[ijkToid(i, j - 1, k)] <= 0 ||
          phi[id] * phi[ijkToid(i, j + 1, k)] <= 0 ||
          phi[id] * phi[ijkToid(i, j, k - 1)] <= 0 ||
          phi[id] * phi[ijkToid(i, j, k + 1)] <= 0) {
        _surfaceCells[id] = true;
#pragma omp critical
        { surface.emplace_back(id); }
      }
    }
  }
  return surface;
}

} // namespace Leviathan