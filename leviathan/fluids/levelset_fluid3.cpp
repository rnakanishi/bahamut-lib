#include <fluids/levelset_fluid3.h>
#include <blas/interpolator.h>
#include <blas/weno.h>
#include <cmath>
#include <iostream>
#include <set>
#include <queue>

namespace Leviathan {

LevelSetFluid3::LevelSetFluid3()
    : LevelSetFluid3(Eigen::Array3i(32, 32, 32),
                     Ramuh::BoundingBox3(Eigen::Array3d(-1, -1, -1),
                                         Eigen::Array3d(1, 1, 1))) {}

LevelSetFluid3::LevelSetFluid3(Eigen::Array3i gridSize,
                               Ramuh::BoundingBox3 domain)
    : Ramuh::MacGrid3(domain, gridSize) {
  _phiId = newCellScalarLabel("phi", 1e8);
  _cellVelocityId = newCellArrayLabel("cellVelocity");
  _faceVelocityId = newFaceScalarLabel("faceVelocity");
  //  --> Even though velocities are vector, MAC grid split them into scalar
  //  pieces for each coordinate
  _originalDt = _dt = 1. / 60;
  _ellapsedDt = 0;
  _tolerance = 1e-10;

  _cellGradientId = newCellArrayLabel("cellGradient");
  _isSurfaceCell.resize(cellCount(), false);
}

void LevelSetFluid3::computeCellsGradient() {
  auto &phi = getCellScalarData("phi");
  auto &gradient = getCellArrayData("cellGradient");
  auto &cellVelocity = getCellArrayData(_cellVelocityId);
  auto h = getH();

#pragma omp parallel for
  for (size_t cellId = 0; cellId < cellCount(); cellId++) {
    auto ijk = idToijk(cellId);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    for (size_t d = 0; d < 3; d++) {
      double velocity = cellVelocity[cellId][d];
      int neighId;

      if (velocity <= 0) {
        if (d == 0)
          neighId = ijkToid(std::max(0, i - 1), j, k);
        else if (d == 1)
          neighId = ijkToid(i, std::max(0, j - 1), k);
        else if (d == 2)
          neighId = ijkToid(i, j, std::max(0, k - 1));

        gradient[cellId][d] = phi[cellId] - phi[neighId];
      } else {
        if (d == 0)
          neighId = ijkToid(std::min(_gridSize[d], i + 1), j, k);
        else if (d == 1)
          neighId = ijkToid(i, std::min(_gridSize[d], j + 1), k);
        else if (d == 2)
          neighId = ijkToid(i, j, std::min(_gridSize[d], k + 1));

        gradient[cellId][d] = phi[neighId] - phi[cellId];
      }
    }
  }
}

void LevelSetFluid3::advectSemiLagrangean() {
  size_t nCells = cellCount();
  auto &phi = getCellScalarData(_phiId);
  auto &cellVelocity = getCellArrayData(_cellVelocityId);
  std::vector<double> newPhi;
  newPhi.insert(newPhi.begin(), phi.begin(), phi.end());

#pragma omp parallel for
  for (int sId = 0; sId < _surfaceCellIds.size(); sId++) {
    int id = _surfaceCellIds[sId];
    Eigen::Array3d position, h, velocity;
    h = getH();
    velocity = cellVelocity[id];

    position = getCellPosition(id);
    position = position - (velocity * _dt);
    Eigen::Array3i index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
    // Check if inside domain
    if ((velocity.matrix().norm() > _tolerance) && index[0] >= 0 &&
        index[1] >= 0 && index[2] >= 0 && index[0] < _gridSize[0] &&
        index[1] < _gridSize[1] && index[2] < _gridSize[2]) {
      newPhi[id] = interpolateCellScalarData(_phiId, position);
    }
  }

#pragma omp parallel for
  for (int sId = 0; sId < _surfaceCellIds.size(); sId++) {
    int id = _surfaceCellIds[sId];
    phi[id] = newPhi[id];
  }
}

void LevelSetFluid3::advectUpwind() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  auto &cellVelocity = getCellArrayData(_cellVelocityId);
  auto &cellGradient = getCellArrayData(_cellGradientId);

#pragma omp parallel for
  // for (int id = 0; id < cellCount(); id++) {
  for (int sId = 0; sId < _surfaceCellIds.size(); sId++) {
    int id = _surfaceCellIds[sId];

    Eigen::Vector3d velocity = cellVelocity[id].matrix();
    Eigen::Vector3d gradient = cellGradient[id].matrix();

    phi[id] = phi[id] - velocity.dot(gradient) * _dt;
  }
}

void LevelSetFluid3::computeWenoGradient() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  auto &gradient = getCellArrayData(_cellGradientId);
  auto &cellVelocity = getCellArrayData(_cellVelocityId);

  computeCellVelocity();
// Weno computation
#pragma omp parallel for
  // for (int id = 0; id < cellCount(); id++) {
  for (int sId = 0; sId < _surfaceCellIds.size(); sId++) {
    int id = _surfaceCellIds[sId];
    std::vector<double> values(6);
    auto ijk = idToijk(id);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    gradient[id] = 0;

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
      velocity = cellVelocity[id][coord];

      // For positive velocity  upwind-based scheme is used. Otherwise,
      // downwind
      if (velocity <= 0) {
        // ========= DOWNWIND =========
        for (int ival = 0, inc = 3; ival < 7; ival++, inc--) {
          int index, index1, faceIndex;
          // Each coordinate have to be treated separatedly due to cell-face
          // indexation
          switch (coord) {
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
      gradient[id][coord] = dPhi;
    }
  }
}

void LevelSetFluid3::advectWeno() {
  auto &gradient = getCellArrayData(_cellGradientId);
  auto &cellVelocity = getCellArrayData(_cellVelocityId);
  auto &phi = getCellScalarData(_phiId);

  findSurfaceCells(8);

  std::vector<double> phiN(phi.begin(), phi.end());
#pragma omp parallel
  {
    computeCellVelocity();
    computeWenoGradient();
    advectUpwind();

    // computeCellVelocity();
    computeWenoGradient();

    advectUpwind();
#pragma omp for
    for (int sId = 0; sId < _surfaceCellIds.size(); sId++) {
      int i = _surfaceCellIds[sId];
      phi[i] = 0.75 * phiN[i] + 0.25 * phi[i];
    }

    // computeCellVelocity();
    computeWenoGradient();
    advectUpwind();

#pragma omp for
    for (int sId = 0; sId < _surfaceCellIds.size(); sId++) {
      int i = _surfaceCellIds[sId];
      phi[i] = phiN[i] / 3 + 2 * phi[i] / 3;
    }
  }
}

void LevelSetFluid3::computeCellVelocity() {
  auto &cellVelocity = getCellArrayData(_cellVelocityId);
  auto h = getH();

#pragma omp parallel for
  for (int sId = 0; sId < _surfaceCellIds.size(); sId++) {
    int id = _surfaceCellIds[sId];
    auto position = getCellPosition(id);

    for (int face = 0; face < 3; face++) {
      auto &data = getFaceScalarData(face, _faceVelocityId);
      auto cellId = idToijk(id);
      int index;
      switch (face) {
      case 0:
        index = faceijkToid(face, cellId[0] + 1, cellId[1], cellId[2]);
        break;
      case 1:
        index = faceijkToid(face, cellId[0], cellId[1] + 1, cellId[2]);
        break;
      case 2:
        index = faceijkToid(face, cellId[0], cellId[1], cellId[2] + 1);
        break;
      default:
        break;
      }
      cellVelocity[id][face] =
          0.5 * (data[faceijkToid(face, cellId[0], cellId[1], cellId[2])] +
                 data[index]);
    }
  }
}

void LevelSetFluid3::advectCip() {
  // FIXME: LevelSetFluid3::advectCip Not working
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  auto &G = getCellArrayData(_cellGradientId);
  auto &uVelocity = getFaceScalarData(0, _faceVelocityId);
  auto &vVelocity = getFaceScalarData(1, _faceVelocityId);
  auto &wVelocity = getFaceScalarData(2, _faceVelocityId);
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
  std::vector<double> gradient(cellCount());
  std::vector<double> cellSignal(cellCount());
  std::vector<double> initialPhi(cellCount());
  std::vector<double> interfaceFactor(cellCount());
  std::vector<bool> isInterface(cellCount(), false);

  double totalError;
  int iterations = 0;
  bool shouldStop = false;
#pragma omp parallel
  {
#pragma omp for schedule(static) nowait
    for (int id = 0; id < cellCount(); id++) {
      initialPhi[id] = phi[id];
      if (phi[id] > 0)
        cellSignal[id] = 1;
      else if (phi[id] < 0)
        cellSignal[id] = -1;
      else
        cellSignal[id] = 0;

      auto ijk = idToijk(id);
      isInterface[id] = false;
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
          interfaceFactor[id] = h[0] * phi[id];

          // These differential values are used to improbe robustness when
          // redistancing levelsets that may suffer too much from topological
          // changes
          double Dphi, dCentral[3], dUp[3], dDown[3];
          dCentral[0] = phi[ijkToid(i + 1, j, k)] - phi[ijkToid(i - 1, j, k)];
          dCentral[1] = phi[ijkToid(i, j + 1, k)] - phi[ijkToid(i, j - 1, k)];
          dCentral[2] = phi[ijkToid(i, j, k + 1)] - phi[ijkToid(i, j, k - 1)];
          dDown[0] = phi[ijkToid(i, j, k)] - phi[ijkToid(i - 1, j, k)];
          dDown[1] = phi[ijkToid(i, j, k)] - phi[ijkToid(i, j - 1, k)];
          dDown[2] = phi[ijkToid(i, j, k)] - phi[ijkToid(i, j, k - 1)];
          dUp[0] = phi[ijkToid(i + 1, j, k)] - phi[ijkToid(i, j, k)];
          dUp[1] = phi[ijkToid(i, j + 1, k)] - phi[ijkToid(i, j, k)];
          dUp[2] = phi[ijkToid(i, j, k + 1)] - phi[ijkToid(i, j, k)];

          for (size_t d = 0; d < 3; d++) {
            dCentral[d] = dCentral[d] * dCentral[d];
            dUp[d] = dUp[d] * dUp[d];
            dDown[d] = dDown[d] * dDown[d];
          }

          Dphi = std::max(
              h[0], std::max(sqrt(dCentral[0] + dCentral[1] + dCentral[2]) / 2,
                             std::max(sqrt(dUp[0] + dUp[1] + dUp[2]),
                                      sqrt(dDown[0] + dDown[1] + dDown[2]))));
          interfaceFactor[id] /= Dphi;
        }
      }
    }
    for (int t = 0; t < 25 && !shouldStop; t++) {
      double error = 0.;
#pragma omp single nowait
      {
        totalError = 0;
        iterations++;
      }
#pragma omp for schedule(static)
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
        // // Time integration step. Euler method
        double newPhi = 0.;

        if (!isInterface[id]) {
          newPhi = phi[id] - dt * cellSignal[id] * gradient[id];
        } else {
          newPhi = phi[id] - (dt / h[0]) * (cellSignal[id] * abs(phi[id]) -
                                            interfaceFactor[id]);
        }
        double newError = abs(phi[id] - newPhi);
        error += newError;
        phi[id] = newPhi;
      }
#pragma omp atomic
      totalError += error;

#pragma omp single
      if (totalError / cellCount() < dt * h[0] * h[0]) {
        shouldStop = true;
      }
    }
  }
  std::cerr << "Run redistance for " << iterations << " iterations\n";
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

  std::fill(_isSurfaceCell.begin(), _isSurfaceCell.end(), false);
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
        _isSurfaceCell[id] = true;
#pragma omp critical
        { surface.emplace_back(id); }
      }
    }
  }
  return surface;
}

std::vector<int> LevelSetFluid3::findSurfaceCells() {
  return findSurfaceCells(1);
}

std::vector<int> LevelSetFluid3::findSurfaceCells(int surfaceDistance) {
  std::vector<int> distanceToSurface(cellCount(), 1e8);
  std::vector<int> visited(cellCount(), false);
  std::set<int> toSeed;

  std::queue<int> cellQueue;
  auto surfaceCells = trackSurface();
  for (auto cell : surfaceCells) {
    distanceToSurface[cell] = 0;
    cellQueue.push(cell);
  }

  // For every cell surface tracked before, compute bfs and mark those cells
  // that are at most 3 cells away from surface
  while (!cellQueue.empty()) {
    int cell = cellQueue.front();
    toSeed.insert(cell);

    cellQueue.pop();
    if (!visited[cell]) {
      visited[cell] = true;
      auto ijk = idToijk(cell);
      int i = ijk[0], j = ijk[1], k = ijk[2];
      int distance = distanceToSurface[cell];

      // All neighbors are taken into consideration
      int neighbors[6];
      neighbors[0] = ijkToid(std::max(0, i - 1), j, k);
      neighbors[1] =
          ijkToid(std::min(LevelSetFluid3::_gridSize[0] - 1, i + 1), j, k);
      neighbors[2] = ijkToid(i, std::max(0, j - 1), k);
      neighbors[3] =
          ijkToid(i, std::min(LevelSetFluid3::_gridSize[1] - 1, j + 1), k);
      neighbors[4] = ijkToid(i, j, std::max(0, k - 1));
      neighbors[5] =
          ijkToid(i, j, std::min(LevelSetFluid3::_gridSize[2] - 1, k + 1));
      for (auto neighbor : neighbors) {
        distanceToSurface[neighbor] =
            std::min(distanceToSurface[neighbor], distance + 1);
        if (!visited[neighbor]) {
          if (distanceToSurface[neighbor] <= surfaceDistance)
            cellQueue.push(neighbor);
        }
      }
    }
  }
  _surfaceCellIds.clear();
  _surfaceCellIds.insert(_surfaceCellIds.begin(), toSeed.begin(), toSeed.end());
  return _surfaceCellIds;
}

} // namespace Leviathan