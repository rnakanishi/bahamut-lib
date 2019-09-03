#include <fluids/levelset_fluid2.h>
#include <fstream>
#include <sstream>
#include <iostream>

namespace Leviathan {

LevelSetFluid2::LevelSetFluid2()
    : LevelSetFluid2(Eigen::Array2i(32, 32), Ramuh::BoundingBox2::unitBox()) {}
LevelSetFluid2::LevelSetFluid2(Eigen::Array2i gridSize,
                               Ramuh::BoundingBox2 domain)
    : MacGrid2(domain, gridSize) {
  _phiId = newScalarLabel("phi", 1e8);
  _velocityId = newFaceScalarLabel("velocity");
  //  --> Even though velocities are vector, MAC grid split them into scalar
  //  pieces for each coordinate
  _dt = 1. / 60;
  _tolerance = 1e-10;

  newArrayLabel("cellGradient");
}

void LevelSetFluid2::computeCellsGradient() {
  auto &phi = getScalarData("phi");
  auto &gradient = getArrayData("cellGradient");
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

void LevelSetFluid2::print() {}

void LevelSetFluid2::redistance() {
  auto h = getH();
  auto &phi = getScalarData(_phiId);
  double eps = h[0], error = 1e8, dt = 0.5 * h[0];
  std::vector<double> gradient(cellCount());
  std::vector<double> cellSignal(cellCount());
  std::vector<double> initialPhi(cellCount());
  std::vector<double> interfaceFactor(cellCount());
  std::vector<bool> isInterface(cellCount(), false);

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
        interfaceFactor[id] = h[0] * 2. * phi[id];
        double dx, dy;
        dx = phi[ijToid(i + 1, j)] - phi[ijToid(i - 1, j)];
        dy = phi[ijToid(i, j + 1)] - phi[ijToid(i, j - 1)];
        interfaceFactor[id] /= sqrt(dx * dx + dy * dy);
        // cellSignal[id] = interfaceFactor[id] / h[0];
      }
    }
  }
  int t;
  for (t = 0; t < 300; t++) {
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
    print();
    if (error / cellCount() < dt * h[0] * h[0])
      break;
  }
}

} // namespace Leviathan
