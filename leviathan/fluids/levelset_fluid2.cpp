#include <fluids/levelset_fluid2.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <blas/weno.h>
#include <blas/interpolator.h>

namespace Leviathan {

LevelSetFluid2::LevelSetFluid2()
    : LevelSetFluid2(Eigen::Array2i(32, 32), Ramuh::BoundingBox2::unitBox()) {}

LevelSetFluid2::LevelSetFluid2(Eigen::Array2i gridSize,
                               Ramuh::BoundingBox2 domain)
    : MacGrid2(domain, gridSize) {
  _phiId = newCellScalarLabel("phi", 1e8);
  _velocityId = newFaceScalarLabel("velocity");
  //  --> Even though velocities are vector, MAC grid split them into scalar
  //  pieces for each coordinate
  _originalDt = _dt = 1. / 60;
  _ellapsedDt = 0.;

  _tolerance = 1e-10;

  _gradientId = newCellArrayLabel("cellGradient");
  _cellVelocityId = newCellArrayLabel("cellVelocity");

  _surfaceCells.resize(cellCount(), false);
}

void LevelSetFluid2::computeCellsGradient() {
  auto &phi = getCellScalarData("phi");
  auto &gradient = getCellArrayData("cellGradient");
  auto h = getH();

#pragma omp parallel for
  for (int j = 0; j < _gridSize[1]; j++)
    for (int i = 0; i < _gridSize[0]; i++) {
      int id = ijToid(i, j);
      if (i < _gridSize[0])
        gradient[id][0] = (phi[ijToid(i + 1, j)] - phi[id]) / h[0];
      else
        gradient[id][0] = (phi[id] - phi[ijToid(i - 1, j)]) / h[0];

      if (j < _gridSize[1])
        gradient[id][1] = (phi[ijToid(i, j + 1)] - phi[id]) / h[1];
      else
        gradient[id][1] = (phi[id] - phi[ijToid(i, j - 1)]) / h[1];
    }
}

void LevelSetFluid2::advectRungeKutta3() {
  auto &phi = getCellScalarData(_phiId);

  std::vector<double> lastPhi(phi.size());
  for (size_t i = 0; i < phi.size(); i++) {
    lastPhi[i] = phi[i];
  }
  advectUpwind();
  advectUpwind();
#pragma omp parallel for
  for (size_t i = 0; i < phi.size(); i++) {
    phi[i] = 0.75 * lastPhi[i] + 0.25 * phi[i];
  }
  advectUpwind();
#pragma omp parallel for
  for (size_t i = 0; i < phi.size(); i++) {
    phi[i] = lastPhi[i] / 3 + 2 * phi[i] / 3;
  }
}

void LevelSetFluid2::advectUpwind() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  std::vector<double> newPhi(cellCount());

  for (int id = 0; id < cellCount(); id++) {
    auto p = getCellPosition(id);
    auto ij = idToij(id);
    int i = ij[0], j = ij[1];
    newPhi[id] = 0;
    for (size_t coord = 0; coord < 2; coord++) {
      auto &uv = getFaceScalarData(coord, _velocityId);
      auto faceij = ij;
      auto faceid = faceijToid(coord, i, j);
      int facei = faceij[0], facej = faceij[1];
      double velocity;
      int cellij = ij[coord];

      switch (coord) {
      case 0:
        velocity = (uv[faceijToid(coord, facei, facej)] +
                    uv[faceijToid(coord, facei + 1, facej)]) /
                   2;
        break;
      case 1:
        velocity = (uv[faceijToid(coord, facei, facej)] +
                    uv[faceijToid(coord, facei, facej + 1)]) /
                   2;
        break;
      }
      int index;
      if (velocity <= 0) {
        switch (coord) {
        case 0:
          index = std::min(_gridSize[coord] - 1, std::max(0, i - 1));
          newPhi[id] +=
              velocity * (phi[ijToid(i, j)] - phi[ijToid(index, j)]) / h[coord];
          break;
        case 1:
          index = std::min(_gridSize[coord] - 1, std::max(0, j - 1));
          newPhi[id] +=
              velocity * (phi[ijToid(i, j)] - phi[ijToid(i, index)]) / h[coord];
          break;
        }
      } else {
        switch (coord) {
        case 0:
          index = std::min(_gridSize[coord] - 1, std::max(1, i + 1));
          newPhi[id] +=
              velocity * (phi[ijToid(index, j)] - phi[ijToid(i, j)]) / h[coord];
          break;
        case 1:
          index = std::min(_gridSize[coord] - 1, std::max(1, j + 1));
          newPhi[id] +=
              velocity * (phi[ijToid(i, index)] - phi[ijToid(i, j)]) / h[coord];
          break;
        }
      }
    }
  }
  for (int id = 0; id < cellCount(); id++) {
    phi[id] = phi[id] - newPhi[id] * _dt;
  }
}

void LevelSetFluid2::advectWeno() {
  auto &gradient = getCellArrayData(_gradientId);
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

void LevelSetFluid2::computeCentralGradient() {
  auto &gradient = getCellArrayData(_gradientId);
  auto &cellVelocity = getCellArrayData(_cellVelocityId);
  auto &phi = getCellScalarData(_phiId);
  auto h = getH();

  computeCellVelocity();
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    auto ij = idToij(id);
    int i = ij[0], j = ij[1];
    Eigen::Array2d velocity = cellVelocity[id];

    for (size_t coord = 0; coord < 2; coord++) {
      switch (coord) {
      case 0:
        int leftId, rightId;
        leftId = std::max(0, i - 1);
        rightId = std::min(_gridSize[0] - 1, i + 1);
        gradient[id][0] = (phi[ijToid(rightId, j)] - phi[ijToid(leftId, j)]) /
                          ((rightId - leftId) * h[0]);
        break;
      case 1:
        int bottomId, topId;
        bottomId = std::max(0, j - 1);
        topId = std::min(_gridSize[1] - 1, j + 1);
        gradient[id][1] = (phi[ijToid(i, topId)] - phi[ijToid(i, bottomId)]) /
                          ((topId - bottomId) * h[1]);
        break;
      }
    }
  }
}

void LevelSetFluid2::advectSemiLagrangian() {
  computeCellVelocity();
  auto &cellVelocity = getCellArrayData(_cellVelocityId);
  auto &phi = getCellScalarData(_phiId);

  std::vector<double> newPhi(cellCount());
  for (size_t cellId = 0; cellId < cellCount(); cellId++) {
    auto position = getCellPosition(cellId);
    position -= cellVelocity[cellId] * _dt;
    newPhi[cellId] = interpolateCellScalarData(_phiId, position);
  }
  for (size_t cellId = 0; cellId < cellCount(); cellId++) {
    phi[cellId] = newPhi[cellId];
  }
}

void LevelSetFluid2::computeWenoGradient() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  auto &gradient = getCellArrayData(_gradientId);
  auto &cellVelocity = getCellArrayData(_cellVelocityId);

  // Weno computation
  for (int id = 0; id < cellCount(); id++) {
    auto p = getCellPosition(id);
    auto ij = idToij(id);
    int i = ij[0], j = ij[1];
    gradient[id] = 0;
    std::vector<double> values(6);

    for (size_t coord = 0; coord < 2; coord++) {
      auto &uv = getFaceScalarData(coord, _velocityId);
      auto faceij = ij;
      auto faceid = faceijToid(coord, i, j);
      int facei = faceij[0], facej = faceij[1];
      double velocity;
      int cellij = ij[coord];

      switch (coord) {
      case 0:
        velocity = (uv[faceijToid(coord, facei, facej)] +
                    uv[faceijToid(coord, facei + 1, facej)]) /
                   2;
        cellVelocity[id][0] = velocity;
        break;
      case 1:
        velocity = (uv[faceijToid(coord, facei, facej)] +
                    uv[faceijToid(coord, facei, facej + 1)]) /
                   2;
        cellVelocity[id][1] = velocity;
        break;
      }

      if (velocity <= 0) {
        for (int ival = 0, inc = 3; ival < 7; ival++, inc--) {
          int index, index1, faceindex;
          switch (coord) {
          case 0:
            index = std::min(_gridSize[coord] - 1, std::max(1, i + inc));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, i + inc - 1));
            faceindex =
                std::min(_gridSize[coord] - 1, std::max(0, facei + inc));
            // velocity = .5 * (uv[faceijToid(coord, faceindex, facej)] +
            //                  uv[faceijToid(coord, faceindex - 1, facej)]);
            values[ival] =
                (phi[ijToid(index, j)] - phi[ijToid(index1, j)]) / h[coord];
            break;
          case 1:
            index = std::min(_gridSize[coord] - 1, std::max(1, j + inc));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, j + inc - 1));
            faceindex =
                std::min(_gridSize[coord] - 1, std::max(0, facej + inc));
            // velocity = .5 * (uv[faceijToid(coord, facei, faceindex)] +
            //                  uv[faceijToid(coord, facei, faceindex - 1)]);
            values[ival] =
                (phi[ijToid(i, index)] - phi[ijToid(i, index1)]) / h[coord];
          default:
            break;
          }
        }
        // Boundary conditions. This ensure weno doesnt use outside values
        if (cellij <= 1)
          values[5] = 1e4;
        if (cellij <= 0)
          values[4] = 1e4;
        if (cellij >= _gridSize[coord] - 2)
          values[0] = 1e4;
        if (cellij >= _gridSize[coord] - 1)
          values[1] = 1e4;
      } else {
        for (int ival = 0, inc = -3; ival < 7; ival++, inc++) {
          int index, index1, faceIndex;
          // Each coordinate have to be treated separatedly due to cell-face
          // indexation
          switch (coord) {
          case 0:
            index = std::min(_gridSize[coord] - 1, std::max(1, i + inc + 1));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, i + inc));
            faceIndex =
                std::min(_gridSize[coord] - 1, std::max(0, facei + inc + 1));
            // velocity = .5 * (uv[faceijToid(coord, faceIndex, facej)] +
            //                  uv[faceijToid(coord, faceIndex - 1, facej)]);
            values[ival] =
                (phi[ijToid(index, j)] - phi[ijToid(index1, j)]) / h[coord];
            break;
          case 1:
            index = std::min(_gridSize[coord] - 1, std::max(1, j + inc + 1));
            index1 = std::min(_gridSize[coord] - 2, std::max(0, j + inc));
            faceIndex =
                std::min(_gridSize[coord] - 1, std::max(0, facej + inc + 1));
            // velocity = .5 * (uv[faceijToid(coord, facei, faceIndex)] +
            //                  uv[faceijToid(coord, facei, faceIndex - 1)]);
            values[ival] =
                (phi[ijToid(i, index)] - phi[ijToid(i, index1)]) / h[coord];
            break;
          }
        }
        if (cellij <= 1)
          values[0] = 1e4;
        if (cellij <= 0)
          values[1] = 1e4;
        if (cellij >= _gridSize[coord] - 2)
          values[5] = 1e4;
        if (cellij >= _gridSize[coord] - 1)
          values[4] = 1e4;
      }
      double dPhi = Ramuh::Weno::evaluate(values, h[coord], false);
      gradient[id][coord] = dPhi;
    }
  }
}

void LevelSetFluid2::advectCip() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  auto &gradient = getCellArrayData(_gradientId);
  auto &uVelocity = getFaceScalarData(0, _velocityId);
  auto &vVelocity = getFaceScalarData(1, _velocityId);
  static int count = 0;
  std::vector<double> a(cellCount()), b(cellCount());
  std::vector<double> newPhi(cellCount(), 0);
  std::vector<Eigen::Array2d> newGrad(cellCount(), Eigen::Array2d(0));

  // TODO: Max velocity for cfl

  for (size_t id = 0; id < cellCount(); id++) {
    auto ij = idToij(id);
    int i = ij[0], j = ij[1];

    double uCellVelocity = 0.5 * (uVelocity[faceijToid(0, i, j)] +
                                  uVelocity[faceijToid(0, i + 1, j)]);
    double vCellVelocity = 0.5 * (vVelocity[faceijToid(1, i, j)] +
                                  vVelocity[faceijToid(1, i, j + 1)]);

    double XX = -uCellVelocity * _dt;
    double YY = -vCellVelocity * _dt;
    int isgn = (uCellVelocity >= 0) ? 1 : -1;
    int jsgn = (vCellVelocity >= 0) ? 1 : -1;
    int im1 = std::max(0, std::min(_gridSize[0] - 1, i - isgn));
    int jm1 = std::max(0, std::min(_gridSize[1] - 1, j - jsgn));
    // double gx[3], gy[3]; // 0: g[i,j]   1: g[i-sgn,j]   2: g[i,j-sgn]

    double A8 = phi[id] - phi[ijToid(im1, j)] - phi[ijToid(i, jm1)] +
                phi[ijToid(im1, jm1)];
    double tmp = gradient[ijToid(im1, j)][1] - gradient[id][1];
    double A1 = ((gradient[ijToid(im1, j)][0] + gradient[id][0]) * h[0] * isgn -
                 2 * (phi[id] - phi[ijToid(im1, j)])) /
                (h[0] * h[0] * h[0] * isgn);
    double A2 =
        (-A8 - (gradient[ijToid(i, jm1)][0] - gradient[id][0]) * h[0] * isgn) /
        (h[0] * h[0] * h[1] * jsgn);
    double A3 =
        (3 * (phi[ijToid(im1, j)] - phi[id]) +
         (gradient[ijToid(im1, j)][0] + 2 * gradient[id][0]) * h[0] * isgn) /
        (h[0] * h[0]);
    double A4 = (A2 * h[0] * h[0] - tmp) / (h[0] * isgn);
    double A5 =
        (-2 * (phi[id] - phi[ijToid(i, jm1)]) +
         (gradient[ijToid(i, jm1)][1] + gradient[id][1]) * h[1] * jsgn) /
        (h[1] * h[1] * h[1] * jsgn);
    double A6 = (-A8 - tmp * h[1] * jsgn) / (h[0] * h[1] * h[1] * isgn);
    double A7 =
        (3 * (phi[ijToid(i, jm1)] - phi[id]) +
         (gradient[ijToid(i, jm1)][1] + 2 * gradient[id][1]) * h[1] * jsgn) /
        (h[1] * h[1]);
    newPhi[id] =
        ((A1 * XX + A2 * YY + A3) * XX + A4 * YY + gradient[id][0]) * XX +
        ((A5 * YY + A6 * XX + A7) * YY + gradient[id][1]) * YY + phi[id];
    newGrad[id][0] = (3 * A1 * XX + 2 * (A2 * YY + A3)) * XX +
                     (A4 + A6 * YY) * YY + gradient[id][0];
    newGrad[id][1] = (3 * A5 * YY + 2 * (A6 * XX + A7)) * YY +
                     (A4 + A2 * XX) * XX + gradient[id][1];
  }
  for (size_t id = 0; id < cellCount(); id++) {
    phi[id] = newPhi[id];
    gradient[id] = newGrad[id];
  }
  // computeCellsGradient();
}

void LevelSetFluid2::computeCellVelocity() {
  auto &cellVelocity = getCellArrayData(_cellVelocityId);

#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    auto ij = idToij(id);
    int i = ij[0], j = ij[1];

    for (size_t coord = 0; coord < 2; coord++) {
      auto &uv = getFaceScalarData(coord, _velocityId);
      double velocity;
      switch (coord) {
      case 0:
        velocity =
            (uv[faceijToid(coord, i, j)] + uv[faceijToid(coord, i + 1, j)]) / 2;
        cellVelocity[id][0] = velocity;
        break;
      case 1:
        velocity =
            (uv[faceijToid(coord, i, j)] + uv[faceijToid(coord, i, j + 1)]) / 2;
        cellVelocity[id][1] = velocity;
        break;
      }
    }
  }
}

void LevelSetFluid2::print() {}

void LevelSetFluid2::redistance() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  double eps = h[0], error = 1e8, dt = 0.5 * h[0];
  std::vector<double> gradient(cellCount());
  std::vector<double> cellSignal(cellCount());
  std::vector<double> initialPhi(cellCount());
  std::vector<double> interfaceFactor(cellCount());
  std::vector<bool> isInterface(cellCount(), false);

#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    initialPhi[id] = phi[id];

    if (phi[id] > 0)
      cellSignal[id] = 1;
    else if (phi[id] < 0)
      cellSignal[id] = -1;
    else
      cellSignal[id] = 0;

    auto ij = idToij(id);
    isInterface[id] = false;
    int i = ij[0], j = ij[1];
    interfaceFactor[id] = 0;
    if (i > 0 && i < _gridSize[0] - 1 && j > 0 && j < _gridSize[1] - 1) {
      if (phi[id] * phi[ijToid(i - 1, j)] <= 0 ||
          phi[id] * phi[ijToid(i + 1, j)] <= 0 ||
          phi[id] * phi[ijToid(i, j - 1)] <= 0 ||
          phi[id] * phi[ijToid(i, j + 1)] <= 0) {
        isInterface[id] = true;
        interfaceFactor[id] = h[0] * phi[id];
        double Dphi;
        double dxCentral = phi[ijToid(i + 1, j)] - phi[ijToid(i - 1, j)];
        double dyCentral = phi[ijToid(i, j + 1)] - phi[ijToid(i, j - 1)];
        double dxUp = phi[ijToid(i, j)] - phi[ijToid(i - 1, j)];
        double dxDown = phi[ijToid(i + 1, j)] - phi[ijToid(i, j)];
        double dyUp = phi[ijToid(i, j)] - phi[ijToid(i, j - 1)];
        double dyDown = phi[ijToid(i, j + 1)] - phi[ijToid(i, j)];

        Dphi = std::max(
            h[0],
            std::max(sqrt(dxCentral * dxCentral + dyCentral * dyCentral) / 2,
                     std::max(sqrt(dxUp * dxUp + dyUp * dyUp),
                              sqrt(dxDown * dxDown + dyDown * dyDown))));
        interfaceFactor[id] /= Dphi;
        // cellSignal[id] = interfaceFactor[id] / h[0];
      }
    }
  }
  int t;
  for (t = 0; t < 300; t++) {
    bool hasError = false;
    // #pragma omp parallel for
    for (int id = 0; id < cellCount(); id++) {
      auto ij = idToij(id);
      int i = ij[0], j = ij[1];
      double dx[2], dy[2]; // Index 0: Dx-, Index 1: Dx+
      dx[0] = dx[1] = dy[0] = dy[1] = 0.;

      if (i > 0)
        dx[0] = (phi[id] - phi[ijToid(i - 1, j)]) / h[0];
      if (i < _gridSize[0] - 1)
        dx[1] = (phi[ijToid(i + 1, j)] - phi[id]) / h[0];

      if (j > 0)
        dy[0] = (phi[id] - phi[ijToid(i, j - 1)]) / h[1];
      if (j < _gridSize[1] - 1)
        dy[1] = (phi[ijToid(i, j + 1)] - phi[id]) / h[1];

      double gx, gy;
      if (initialPhi[id] > 0) {
        double a, b;
        a = std::max(0., dx[0]);
        b = std::min(0., dx[1]);
        gx = std::max(a * a, b * b);
        a = std::max(0., dy[0]);
        b = std::min(0., dy[1]);
        gy = std::max(a * a, b * b);
      } else if (initialPhi[id] < 0) {
        double a, b;
        a = std::min(0., dx[0]);
        b = std::max(0., dx[1]);
        gx = std::max(a * a, b * b);
        a = std::min(0., dy[0]);
        b = std::max(0., dy[1]);
        gy = std::max(a * a, b * b);
      } else {
        gx = gy = 0.5;
      }
      gradient[id] = std::sqrt(gx + gy) - 1;
    }

    error = 0.;
    // #pragma omp parallel for reduction(+ : error)
    for (int id = 0; id < cellCount(); id++) {
      auto ij = idToij(id);
      int i = ij[0], j = ij[1];
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
}

void LevelSetFluid2::redistanceSimple() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  double eps = h[0], error = 1e8, dt = 0.5 * h[0];
  auto &cellGradient = getCellArrayData(_gradientId);
  std::vector<double> gradient(cellCount());
  std::vector<double> cellSignal(cellCount());
  std::vector<double> initialPhi(cellCount());
  for (int id = 0; id < cellCount(); id++) {
    initialPhi[id] = phi[id];
    cellSignal[id] = phi[id] / sqrt(phi[id] * phi[id] + eps * eps);
  }
  int t;
  for (t = 0; t < 300; t++) {
    bool hasError = false;
    // #pragma omp parallel for
    for (int id = 0; id < cellCount(); id++) {
      auto ij = idToij(id);
      int i = ij[0], j = ij[1];
      double dx[2], dy[2]; // Index 0: Dx-, Index 1: Dx+
      dx[0] = dx[1] = dy[0] = dy[1] = 0.;

      if (i > 0)
        dx[0] = (phi[id] - phi[ijToid(i - 1, j)]) / h[0];
      if (i < _gridSize[0] - 1)
        dx[1] = (phi[ijToid(i + 1, j)] - phi[id]) / h[0];

      if (j > 0)
        dy[0] = (phi[id] - phi[ijToid(i, j - 1)]) / h[1];
      if (j < _gridSize[1] - 1)
        dy[1] = (phi[ijToid(i, j + 1)] - phi[id]) / h[1];

      double gx, gy;
      if (initialPhi[id] > 0) {
        double a, b;
        a = std::max(0., dx[0]);
        b = std::min(0., dx[1]);
        gx = std::max(a * a, b * b);
        a = std::max(0., dy[0]);
        b = std::min(0., dy[1]);
        gy = std::max(a * a, b * b);
      } else if (initialPhi[id] < 0) {
        double a, b;
        a = std::min(0., dx[0]);
        b = std::max(0., dx[1]);
        gx = std::max(a * a, b * b);
        a = std::min(0., dy[0]);
        b = std::max(0., dy[1]);
        gy = std::max(a * a, b * b);
      } else {
        gx = gy = 0.5;
      }
      gradient[id] = std::sqrt(gx + gy) - 1;
    }

    error = 0.;
    // #pragma omp parallel for reduction(+ : error)
    for (int id = 0; id < cellCount(); id++) {
      auto ij = idToij(id);
      int i = ij[0], j = ij[1];
      double newPhi;
      newPhi = phi[id] - dt * cellSignal[id] * gradient[id];

      error += abs(phi[id] - newPhi);
      if (abs(phi[id] - newPhi) > dt * h[0] * h[0])
        hasError = true;

      phi[id] = newPhi;
    }

    if (error / cellCount() < dt * h[0] * h[0])
      // if (!hasError)
      break;
  }
}

bool LevelSetFluid2::advanceTime() {
  _ellapsedDt += _dt;
  if (_ellapsedDt < _originalDt)
    return false;
  _ellapsedDt = 0.0;
  _dt = _originalDt;
  return true;
}

void LevelSetFluid2::applyCfl() {
  Eigen::Vector2d maxVel = Eigen::Vector2d(0.0, 0.0);
  auto &u = getFaceScalarData(0, _velocityId);
  auto &v = getFaceScalarData(1, _velocityId);
  auto h = getH();

  for (int id = 0; id < cellCount(); id++) {
    auto ijk = idToij(id);
    int i, j;
    i = ijk[0];
    j = ijk[1];

    Eigen::Vector2d vel;
    vel[0] = (u[id] + u[faceijToid(0, i + 1, j)]) / 2.0;
    vel[1] = (v[id] + v[faceijToid(1, i, j + 1)]) / 2.0;

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

std::vector<int> LevelSetFluid2::trackSurface() {
  auto h = getH();
  auto &phi = getCellScalarData(_phiId);
  double eps = h[0], error = 1e8, dt = 0.5 * h[0];
  std::vector<double> cellSignal(cellCount());
  std::vector<bool> isInterface(cellCount(), false);
  std::vector<int> surface;

  std::fill(_surfaceCells.begin(), _surfaceCells.end(), false);
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    auto ij = idToij(id);
    int i = ij[0], j = ij[1];
    if (i > 0 && i < _gridSize[0] - 1 && j > 0 && j < _gridSize[1] - 1) {
      if (phi[id] * phi[ijToid(i - 1, j)] <= 0 ||
          phi[id] * phi[ijToid(i + 1, j)] <= 0 ||
          phi[id] * phi[ijToid(i, j - 1)] <= 0 ||
          phi[id] * phi[ijToid(i, j + 1)] <= 0) {
        _surfaceCells[id] = true;
#pragma omp critical
        { surface.emplace_back(id); }
      }
    }
  }
  return surface;
}

} // namespace Leviathan
