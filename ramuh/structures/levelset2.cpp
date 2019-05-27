#include <structures/levelset2.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <utils/macros.h>
#include <omp.h>
#include <cmath>
#include <queue>
#include <set>
#include <blas/interpolator.h>
#include <fstream>

namespace Ramuh {

LevelSet2::LevelSet2() : RegularGrid2() {
  _phi.resize(_resolution[0]);
  for (auto &row : _phi) {
    row.resize(_resolution[1], 1e8);
  }
  _gradPhi.resize(_resolution[0] + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution[1] + 1);
  }
}

void LevelSet2::setResolution(Eigen::Array2i newResolution) {
  RegularGrid2::setResolution(newResolution);
  _phi.resize(_resolution[0]);
  for (auto &row : _phi) {
    row.resize(_resolution[1], 1e8);
  }
  _gradPhi.resize(_resolution[0] + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution[1] + 1);
  }
  _velocity.resize(_resolution[0] + 1);
  for (auto &row : _velocity) {
    row.resize(_resolution[1] + 1);
  }
}

void LevelSet2::setPressureOrder(bool order) { _isPressure2nd = order; }

void LevelSet2::printVertexVelocity() {
  std::cerr << "==== x component: \n";
  for (int j = 0; j < _resolution[1] + 1; j++) {
    for (int i = 0; i < _resolution[0] + 1; i++) {
      std::cerr << _velocity[i][j].x() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "==== y component: \n";
  for (int j = 0; j < _resolution[1] + 1; j++) {
    for (int i = 0; i < _resolution[0] + 1; i++) {
      std::cerr << _velocity[i][j].y() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
}

void LevelSet2::addSphereSurface(Vector2d center, double radius) {
  addSphereSurface(Eigen::Array2d(center[0], center[1]), radius);
}

void LevelSet2::addSphereSurface(Eigen::Array2d center, double radius) {
  for (int j = 0; j < _resolution.y(); j++)
    for (int i = 0; i < _resolution.x(); i++) {
      Eigen::Array2d h(_h[0], _h[1]);
      Eigen::Array2d position = Eigen::Array2d(i, j).cwiseProduct(h) + h / 2.0;
      double distance = (position - center).matrix().norm() - radius;
      _phi[i][j] = std::min(_phi[i][j], distance);
    }
  checkCellMaterial();
}

void LevelSet2::addCubeSurface(Eigen::Array2d lower, Eigen::Array2d upper) {
  addCubeSurface(Vector2d(lower[0], lower[1]), Vector2d(upper[0], upper[1]));
}

void LevelSet2::addCubeSurface(Vector2d lower, Vector2d upper) {
  Vector2d center = (lower + upper) / 2;
  Vector2d size = upper - lower;
  Vector2d h(_h[0], _h[1]);

  for (int j = 0; j < _resolution.y(); j++)
    for (int i = 0; i < _resolution.x(); i++) {
      Vector2d position(Vector2i(i, j) * h + h / 2);
      Vector2d distanceLower = (position - lower).abs();
      Vector2d distanceUpper = (position - upper).abs();
      double phi;
      if (position > lower && position < upper) {
        // Inside cube
        phi = std::min(std::fabs(_phi[i][j]),
                       std::min(distanceLower.min(), distanceUpper.min()));
        _phi[i][j] = -std::min(std::fabs(_phi[i][j]), phi);
      } else {
        // Compute distance to the cube
        position = position - center;
        position = position * 2 / size;

        position = position.abs();

        phi = 0.0;
        phi += std::max(0.0, position[0] - 1);
        phi += std::max(0.0, position[1] - 1);
        phi = sqrt(phi);
        _phi[i][j] = std::min((_phi[i][j]), phi);
      }
    }
  checkCellMaterial();
}

void LevelSet2::checkCellMaterial() {
  _fluidCells.clear();
  // Mark all cells and faces as air
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array2i ij = idToij(_id);
    int i, j;
    i = ij[0];
    j = ij[1];
    _material[i][j] = Material::FluidMaterial::AIR;
    _uFaceMaterial[i][j] = Material::FluidMaterial::AIR;
    _vFaceMaterial[i][j] = Material::FluidMaterial::AIR;
    _uFaceMaterial[i + 1][j] = Material::FluidMaterial::AIR;
    _vFaceMaterial[i][j + 1] = Material::FluidMaterial::AIR;
  }

  // Walk through all cell/faces again to mark them as fluid
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array2i ij = idToij(_id);
    int i, j;
    i = ij[0];
    j = ij[1];
    if (_phi[i][j] <= 0) {
      _material[i][j] = Material::FluidMaterial::FLUID;
      _uFaceMaterial[i][j] = Material::FluidMaterial::FLUID;
      _vFaceMaterial[i][j] = Material::FluidMaterial::FLUID;
      _uFaceMaterial[i + 1][j] = Material::FluidMaterial::FLUID;
      _vFaceMaterial[i][j + 1] = Material::FluidMaterial::FLUID;

#pragma omp critical
      { _fluidCells.emplace_back(_id); }
    }
  }
}

void LevelSet2::addImplicitFunction() { NOT_IMPLEMENTED(); }

void LevelSet2::macCormackAdvection() {
  std::vector<std::vector<double>> newPhi;
  std::vector<std::vector<bool>> phiChanged;
  newPhi.resize(_resolution.x() + 1);
  for (auto &row : newPhi) {
    row.resize(_resolution.y());
  }
  phiChanged.resize(_resolution.x() + 1);
  for (auto &row : phiChanged) {
    row.resize(_resolution.y());
  }
  Eigen::Array2i resolution(_resolution[0], _resolution[1]);
#pragma omp for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array2i ij = idToij(_id);
    int i, j;
    i = ij[0];
    j = ij[1];

    // if (std::fabs(_phi[i][j]) > 5 * _h[0]) {
    //   newPhi[i][j] = _phi[i][j];
    //   continue;
    // }
    Eigen::Array2d position, h(_h[0], _h[1]);
    Eigen::Vector2d velocity;
    Eigen::Array2i index;
    // Find threshold values for clamping
    double clamp[2]; // 0: min, 1: max
    clamp[0] = 1e8;
    clamp[1] = -1e8;

    // Find the point as if it was in previous timestep and work with it
    position = Eigen::Array2d(i, j).cwiseProduct(h) + (h / 2.0);
    velocity[0] = _interpolateVelocityU(position);
    velocity[1] = _interpolateVelocityV(position);
    position -= velocity.array() * _dt;

    index = Eigen::floor(position.cwiseQuotient(h)).cast<int>();
    // TODO: look for a better solution
    if ((velocity.norm() > _tolerance) && index[0] >= 0 && index[1] >= 0 &&
        index[2] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1] && index[2] < _resolution[2]) {

      double phi_n;
      int phiSignal = (_phi[i][j] >= 0) ? 1 : -1;
      phi_n = _interpolatePhi(position, clamp[0], clamp[1]);
      // Advect forward in time the value there
      velocity[0] = _interpolateVelocityU(position);
      velocity[1] = _interpolateVelocityV(position);
      position += velocity.array() * _dt;
      try {
        double phi_n1_hat = _interpolatePhi(position);
        // double phi_n1_hat = _interpolatePhi(position, phiSignal);
        // Analyse the error from original to the forward advected
        double error = 0.5 * (_phi[i][j] - phi_n1_hat);
        newPhi[i][j] = std::max(clamp[0], std::min(clamp[1], phi_n + error));
        // newPhi[i][j] = phi_n + error;
      } catch (const char *error) {
        std::cerr << error;
        std::cerr << "Velocity " << velocity.transpose() << std::endl;
        std::cerr << "Original phi: " << _phi[i][j] << std::endl;
        std::cerr << "Interpolated phi_n: " << phi_n << std::endl;
        throw("levelSet3::macComarckAdvection\n");
      }
    } else
      newPhi[i][j] = _phi[i][j];
  }
#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      _phi[i][j] = newPhi[i][j];
}

void LevelSet2::integrateLevelSet() {
  // std::vector<std::vector<double>> oldPhi;
  std::vector<std::vector<double>> oldPhi;
  oldPhi.resize(_resolution.x() + 1);
  for (auto &row : oldPhi) {
    row.resize(_resolution.y());
  }

#pragma omp for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array2i ij = idToij(_id);
    int i, j;
    i = ij[0];
    j = ij[1];

    // if (std::fabs(_phi[i][j]) > 5 * _h[0]) {
    //   oldPhi[i][j] = _phi[i][j];
    //   continue;
    // }

    Eigen::Array2d position, backPosition, velocity, h(_h.x(), _h.y()),
        cellCenter;
    Eigen::Array2i index;
    double newPhi = _phi[i][j];

    position = Eigen::Array2d(i, j).cwiseProduct(h) +
               Eigen::Array2d(h[0] / 2, h[1] / 2);

    velocity[0] = _interpolateVelocityU(position);
    velocity[1] = _interpolateVelocityV(position);

    backPosition = position - (velocity * _dt);
    index = Eigen::floor(backPosition.cwiseQuotient(h)).cast<int>();

    // Check if inside domain
    if ((velocity.matrix().norm() > _tolerance) && index[0] >= 0 &&
        index[1] >= 0 && index[0] < _resolution[0] &&
        index[1] < _resolution[1]) {
      newPhi = _interpolatePhi(backPosition);
    }
    oldPhi[i][j] = newPhi;
  }
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array2i ij = idToij(_id);
    int i, j;
    i = ij[0];
    j = ij[1];

    _phi[i][j] = oldPhi[i][j];
  }
}

void LevelSet2::redistance() {
  // TODO: Maybe change this fucntion
  std::vector<std::vector<double>> tempPhi;
  std::vector<std::vector<bool>> processed;
  std::priority_queue<std::pair<double, int>,
                      std::vector<std::pair<double, int>>,
                      std::greater<std::pair<double, int>>>
      cellsQueue;
  std::set<int> cellsAdded;

  tempPhi.resize(_resolution[0]);
  processed.resize(_resolution[0]);
  for (int i = 0; i < _resolution[0]; i++) {
    tempPhi[i].resize(_resolution[1], 1e8);
    processed[i].resize(_resolution[1], false);
  }
#pragma omp parallel for
  // While setting remaining cells to infinity
  for (int i = 0; i < _resolution[0]; i++)
    for (int j = 0; j < _resolution[1]; j++) {
      double cellPhi;
      int cellSign;
      Material::FluidMaterial cellMaterial;
      cellMaterial = _material[i][j];
      cellPhi = _phi[i][j];
      cellSign = cellPhi / std::fabs(cellPhi);

      // Look for surface cells and compute distance over them
      if (i < _resolution[0] - 1 &&
          std::signbit(cellPhi) != std::signbit(_phi[i + 1][j])) {
        // compute distance to the surface
        double theta = cellSign * cellPhi / ((cellPhi) - (_phi[i + 1][j]));
        tempPhi[i][j] = cellSign * std::min(std::fabs(tempPhi[i][j]),
                                            std::fabs(theta * _h.x()));
        processed[i][j] = true;
      }
      if (i > 0 && std::signbit(cellPhi) != std::signbit(_phi[i - 1][j])) {
        double theta = cellSign * cellPhi / ((cellPhi) - (_phi[i - 1][j]));
        tempPhi[i][j] = cellSign * std::min(std::fabs(tempPhi[i][j]),
                                            std::fabs(theta * _h.x()));
        processed[i][j] = true;
      }
      if (j < _resolution[1] - 1 &&
          std::signbit(cellPhi) != std::signbit(_phi[i][j + 1])) {
        double theta = cellSign * cellPhi / ((cellPhi) - (_phi[i][j + 1]));
        tempPhi[i][j] = cellSign * std::min(std::fabs(tempPhi[i][j]),
                                            std::fabs(theta * _h.x()));
        processed[i][j] = true;
      }
      if (j > 0 && std::signbit(cellPhi) != std::signbit(_phi[i][j - 1])) {
        double theta = cellSign * cellPhi / ((cellPhi) - (_phi[i][j - 1]));
        tempPhi[i][j] = cellSign * std::min(std::fabs(tempPhi[i][j]),
                                            std::fabs(theta * _h.x()));
        processed[i][j] = true;
      }

      if (processed[i][j]) {
#pragma omp critical
        {
          cellsQueue.push(
              std::make_pair(std::fabs(tempPhi[i][j]), ijToId(i, j)));
          cellsAdded.insert(ijToId(i, j));
        }
      }
    }

  // Propagating distances
  while (!cellsQueue.empty()) {
    int cellId = cellsQueue.top().second;
    cellsQueue.pop();
    int i = cellId % _resolution[0], j = cellId / _resolution[0];
    // std::cerr << "Id " << cellId << " -> " << i << ' ' << j << std::endl;
    // If the cell is not processed yet, compute its distance
    if (!processed[i][j]) {
      processed[i][j] = true;
      double distances[2] = {1e8, 1e8};
      int distanceSignal = _phi[i][j] / std::fabs(_phi[i][j]);
      distances[0] =
          std::min(std::fabs(tempPhi[std::max(0, i - 1)][j]),
                   std::fabs(tempPhi[std::min(_resolution[0] - 1, i + 1)][j]));
      distances[1] =
          std::min(std::fabs(tempPhi[i][std::max(0, j - 1)]),
                   std::fabs(tempPhi[i][std::min(_resolution[1] - 1, j + 1)]));
      if (distances[0] > distances[1]) {
        double aux = distances[0];
        distances[0] = distances[1];
        distances[1] = aux;
      }
      // std::cerr << "Distances " << distances[0] << ' ' << distances[1]
      // << std::endl;
      double newPhi = distances[0] + _h.x();
      if (std::signbit(distances[0]) == std::signbit(distances[1]) &&
          std::fabs(newPhi) > std::fabs(distances[1])) {
        double h = _h.x() * _h.x();
        newPhi =
            0.5 *
            (distances[0] + distances[1] +
             std::sqrt(2 * _h.x() * _h.x() - ((distances[1] - distances[0]) *
                                              (distances[1] - distances[0]))));
      }
      newPhi *= distanceSignal;
      if (std::fabs(newPhi) < std::fabs(tempPhi[i][j]))
        tempPhi[i][j] = newPhi;
    }
    // Compute neighbor distance to the levelset and add it to queue
    if (i > 0 && cellsAdded.find(ijToId(i - 1, j)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i - 1][j]), ijToId(i - 1, j)));
      cellsAdded.insert(ijToId(i - 1, j));
    }
    if (i < _resolution[0] - 1 &&
        cellsAdded.find(ijToId(i + 1, j)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i + 1][j]), ijToId(i + 1, j)));
      cellsAdded.insert(ijToId(i + 1, j));
    }
    if (j > 0 && cellsAdded.find(ijToId(i, j - 1)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i][j - 1]), ijToId(i, j - 1)));
      cellsAdded.insert(ijToId(i, j - 1));
    }
    if (j < _resolution[1] - 1 &&
        cellsAdded.find(ijToId(i, j + 1)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i][j + 1]), ijToId(i, j + 1)));
      cellsAdded.insert(ijToId(i, j + 1));
    }
  }
  // Solve Eikonal function
  for (int j = 0; j < _resolution[1]; j++) {
    for (int i = 0; i < _resolution[0]; i++) {
      _phi[i][j] = tempPhi[i][j];
    }
  }
}

void LevelSet2::printLevelSetValue() {
  for (int j = _resolution[1] - 1; j >= 0; j--) {
    for (int i = 0; i < _resolution[0]; i++) {
      std::cout << _phi[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

std::vector<double> &LevelSet2::operator[](const int i) { return _phi[i]; }

void LevelSet2::solvePressure() {
  // Compute velocity divergent over cell center
  // std::vector<double> divergent(cellCount(), 0.0);
  // TODO: Change to fluid cells count
  Eigen::VectorXd divergent, pressure;
  Eigen::Array2d h(_h[0], _h[1]);
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
  int ghostId = _getMapId(-1);
#pragma omp parallel for
  for (int _id = 0; _id < cellCount(); _id++) {
    Eigen::Array2i ij = idToij(_id);
    int i, j;
    i = ij[0];
    j = ij[1];

    if (_phi[i][j] > 0.0)
      continue;
    int id = _getMapId(_id);

    int validCells = 0;
    double h2 = _h.x() * _h.x();
    double centerWeight = 0.0;
    std::vector<Eigen::Triplet<double>> threadTriplet;
    threadTriplet.clear();
    Eigen::Array2d center = Eigen::Array2d(i, j).cwiseProduct(h) + h / 2.0;

    if (i > 0) {
      validCells++;
      if (_phi[i - 1][j] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijToId(i - 1, j)), -1 / (h2));
      } else {
        if (_isPressure2nd) {
          double theta;
          // Eigen::Array2d surface = _findSurfaceCoordinate(
          // Eigen::Array2i(i, j), Eigen::Array2i(i - 1, j));
          // theta = (surface - center).matrix().norm() / _h.x();
          theta = -_phi[i][j] / (_phi[i - 1][j] - _phi[i][j]);
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
      if (_phi[i + 1][j] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijToId(i + 1, j)), -1 / (h2));
      } else {
        if (_isPressure2nd) {
          double theta;
          // Eigen::Array2d surface = _findSurfaceCoordinate(
          // Eigen::Array2i(i, j), Eigen::Array2i(i + 1, j));
          // theta = (surface - center).matrix().norm() / _h.x();
          theta = -_phi[i][j] / (_phi[i + 1][j] - _phi[i][j]);
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
      if (_phi[i][j - 1] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijToId(i, j - 1)), -1 / (h2));
      } else {
        if (_isPressure2nd) {
          double theta;
          // Eigen::Array2d surface = _findSurfaceCoordinate(
          // Eigen::Array2i(i, j), Eigen::Array2i(i, j - 1));
          // theta = (surface - center).matrix().norm() / _h.y();
          theta = -_phi[i][j] / (_phi[i][j - 1] - _phi[i][j]);
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
      if (_phi[i][j + 1] < 0.0) {
        centerWeight += 1 / h2;
        threadTriplet.emplace_back(id, _getMapId(ijToId(i, j + 1)), -1 / (h2));
      } else {
        if (_isPressure2nd) {
          double theta;
          // Eigen::Array2d surface = _findSurfaceCoordinate(
          // Eigen::Array2i(i, j), Eigen::Array2i(i, j + 1));
          // theta = (surface - center).matrix().norm() / _h.y();
          theta = -_phi[i][j] / (_phi[i][j + 1] - _phi[i][j]);
          centerWeight += (1 / (h2 * theta));
          threadTriplet.emplace_back(id, ghostId, -1 / (h2 * theta));
        } else {
          centerWeight += (1 / h2);
          threadTriplet.emplace_back(id, ghostId, -1 / (h2));
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
    divergent[id] -= (_u[i + 1][j] - _u[i][j]) / _h[0];
    divergent[id] -= (_v[i][j + 1] - _v[i][j]) / _h[1];
    divergent[id] /= _dt;
  }
  divergent[ghostId] = 0;
  triplets.emplace_back(ghostId, ghostId, 1);

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
  _maxVelocity[0] = _maxVelocity[1] = -1e8;
#pragma omp parallel for
  for (int id = 0; id < cellCount(); id++) {
    Eigen::Array2i ij = idToij(id);
    int i, j;
    i = ij[0];
    j = ij[1];

    if (_material[i][j] != Material::FluidMaterial::FLUID) {
      bool hasFluidNeighbor = false;
      if (i > 0 && _material[i - 1][j] == Material::FluidMaterial::FLUID)
        hasFluidNeighbor = true;
      if (j > 0 && _material[i][j - 1] == Material::FluidMaterial::FLUID)
        hasFluidNeighbor = true;
      if (!hasFluidNeighbor)
        continue;
    }
    double pressure1, pressure2;
    auto it = idMap.find(id);

    pressure1 = (it == idMap.end()) ? 0.0 : pressure[it->second];
    pressure2 = 0.0;
    if (i > 0) {
      auto it = idMap.find(ijToId(i - 1, j));
      if (it == idMap.end())
        pressure2 = 0.0;
      else
        pressure2 = pressure[it->second];
      _u[i][j] -= (_dt * (pressure1 - pressure2) / _h.x());
      if (std::isnan(_u[i][j]) || std::isinf(_u[i][j])) {
        std::cerr << "Infinite velocity component";
      }
      _maxVelocity[0] = std::max(_maxVelocity[0], std::fabs(_u[i][j]));
    }
    pressure2 = 0.0;
    if (j > 0) {
      auto it = idMap.find(ijToId(i, j - 1));
      if (it == idMap.end())
        pressure2 = 0.0;
      else
        pressure2 = pressure[it->second];
      _v[i][j] -= (_dt * (pressure1 - pressure2) / _h.y());
      if (std::isnan(_v[i][j]) || std::isinf(_v[i][j])) {
        std::cerr << "Infinite velocity component";
      }
      _maxVelocity[1] = std::max(_maxVelocity[1], std::fabs(_v[i][j]));
    }
  }
  if (_maxVelocity[0] > 1e4 || _maxVelocity[1] > 1e4) {
    std::cerr << "Vellocity too high\n";
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

double LevelSet2::_interpolatePhi(Eigen::Array2d position, double &min,
                                  double &max) {
  Eigen::Array2i index = Eigen::floor(position.cwiseQuotient(_h)).cast<int>();
  Eigen::Array2d cellCenter = index.cast<double>().cwiseProduct(_h) + _h / 2.0;

  if (position[0] < cellCenter[0])
    index[0]--;
  if (position[1] < cellCenter[1])
    index[1]--;
  std::vector<int> iCandidates, jCandidates;
  for (int i = -1; i < 3; i++) {
    iCandidates.push_back(index[0] + 1);
    jCandidates.push_back(index[1] + 1);
  }
  std::vector<Eigen::Array2d> points;
  std::vector<double> values;
  min = 1e8;
  max = -1e8;

  for (auto v : jCandidates) {
    for (auto u : iCandidates) {
      int x, y;
      points.emplace_back((u + 0.5) * _h[0], (v + 0.5) * _h[1]);
      x = std::max(0, std::min(_resolution[0] - 1, u));
      y = std::max(0, std::min(_resolution[1] - 1, v));
      values.emplace_back(_phi[x][y]);

      min = std::min(values.back(), min);
      max = std::max(values.back(), max);
    }
    double arrayPos[] = {position[0], position[1]};
    double phi = Interpolator::bicubic(arrayPos, points, values);
  }
}

double LevelSet2::_interpolatePhi(Eigen::Array2d position) {
  double min, max;
  _interpolatePhi(position, min, max);
}

TriangleMesh LevelSet2::marchingTriangles() {
  TriangleMesh mesh;

  return mesh;
}

} // namespace Ramuh