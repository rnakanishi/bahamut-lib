#include <structures/levelset2.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <utils/macros.h>
#include <omp.h>
#include <cmath>
#include <queue>
#include <set>

namespace Ramuh {

LevelSet2::LevelSet2() : RegularGrid2() {
  _phi.resize(_resolution.x());
  for (auto &row : _phi) {
    row.resize(_resolution.y(), 1e8);
  }
  _gradPhi.resize(_resolution.x() + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution.y() + 1);
  }
}

void LevelSet2::setResolution(Vector2i newResolution) {
  RegularGrid2::setResolution(newResolution);
  _phi.resize(_resolution.x());
  for (auto &row : _phi) {
    row.resize(_resolution.y(), 1e8);
  }
  _gradPhi.resize(_resolution.x() + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution.y() + 1);
  }
  _velocity.resize(_resolution.x() + 1);
  for (auto &row : _velocity) {
    row.resize(_resolution.y() + 1);
  }
}

void LevelSet2::printVertexVelocity() {
  std::cerr << "==== x component: \n";
  for (int j = 0; j < _resolution.y() + 1; j++) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _velocity[i][j].x() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "==== y component: \n";
  for (int j = 0; j < _resolution.y() + 1; j++) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _velocity[i][j].y() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
}

void LevelSet2::addSphereSurface(Vector2d center, double radius) {
  for (int j = 0; j < _resolution.y(); j++)
    for (int i = 0; i < _resolution.x(); i++) {
      Vector2d position(Vector2i(i, j) * _h + _h / 2);
      double distance = (position - center).length() - radius;
      _phi[i][j] = distance;
    }
  checkCellMaterial();
}

void LevelSet2::addCubeSurface(Vector2d lower, Vector2d upper) {
  for (int j = 0; j < _resolution.y(); j++)
    for (int i = 0; i < _resolution.x(); i++) {
      Vector2d position(Vector2i(i, j) * _h + _h / 2);
      Vector2d distanceLower = (position - lower).abs();
      Vector2d distanceUpper = (position - upper).abs();
      if (position > lower && position < upper) {
        _phi[i][j] = -std::min(
            _phi[i][j], std::min(distanceLower.min(), distanceUpper.min()));
      } else {
        _phi[i][j] = std::min(
            _phi[i][j], std::min(distanceLower.min(), distanceUpper.min()));
      }
    }
  checkCellMaterial();
}

void LevelSet2::checkCellMaterial() {
  for (int j = 0; j < _resolution.y(); j++)
    for (int i = 0; i < _resolution.x(); i++) {
      if (_phi[i][j] <= 0)
        _material[i][j] = Material::FluidMaterial::FLUID;
      else
        _material[i][j] = Material::FluidMaterial::AIR;
    }
}

void LevelSet2::addImplicitFunction() { NOT_IMPLEMENTED(); }

void LevelSet2::integrateLevelSet() {
  int cellCount;
  cellCount = _resolution.x() * _resolution.y();
  // std::vector<std::vector<double>> oldPhi;
  double **oldPhi;
  oldPhi = new double *[_resolution.x()];
  for (int i = 0; i < _resolution.x(); i++) {
    oldPhi[i] = new double[_resolution.y()];
  }

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      oldPhi[i][j] = _phi[i][j];

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++) {
      Vector2d position, backPosition, velocity, h(_h.x(), _h.y()), cellCenter;
      Vector2i index;
      double newPhi = oldPhi[i][j];
      double distanceCount = 0.0, distance = 0.0;

      cellCenter = position = Vector2d(i, j) * h + h / 2.0;
      velocity.x((_u[i][j].x() + _u[i + 1][j].x()) / 2.0);
      velocity.y((_v[i][j].y() + _v[i][j + 1].y()) / 2.0);
      backPosition = position - (velocity * _dt);

      index.x(std::floor(backPosition.x() / h.x()));
      index.y(std::floor(backPosition.y() / h.y()));

      if (velocity.length() > 1e-8 && index >= Vector2i(0, 0) &&
          index < Vector2(_resolution.x(), _resolution.y())) {
        // If inside simulation domain
        cellCenter = Vector2d(index.x(), index.y()) * h + h / 2.0;
        std::vector<int> iCandidates, jCandidates;
        iCandidates.push_back(index.x());
        jCandidates.push_back(index.y());

        if (backPosition.x() > cellCenter.x() &&
            index.x() < _resolution.x() - 1)
          iCandidates.push_back(index.x() + 1);
        else if (backPosition.x() < cellCenter.x() && index.x() > 0)
          iCandidates.push_back(index.x() - 1);

        if (backPosition.y() > cellCenter.y() &&
            index.y() < _resolution.y() - 1)
          jCandidates.push_back(index.y() + 1);
        else if (backPosition.y() < cellCenter.y() && index.y() > 0)
          jCandidates.push_back(index.y() - 1);

        newPhi = 0.;
        for (auto u : iCandidates)
          for (auto v : jCandidates) {
            position = Vector2i(u, v) * h + h / 2;
            distance = (backPosition - position).length();
            distanceCount += 1. / distance;
            newPhi += (1. / distance) * oldPhi[u][v];
          }
        newPhi /= distanceCount;
      }
      _phi[i][j] = newPhi;
    }

  // Clean up
  for (int i = 0; i < _resolution.y(); ++i) {
    delete[] oldPhi[i];
  }
  delete[] oldPhi;
}

void LevelSet2::redistance() {
  std::vector<std::vector<double>> tempPhi;
  std::vector<std::vector<bool>> processed;
  std::priority_queue<std::pair<double, int>,
                      std::vector<std::pair<double, int>>,
                      std::greater<std::pair<double, int>>>
      cellsQueue;
  std::set<int> cellsAdded;

  tempPhi.resize(_resolution.x());
  processed.resize(_resolution.x());
  for (int i = 0; i < _resolution.x(); i++) {
    tempPhi[i].resize(_resolution.y(), 1e8);
    processed[i].resize(_resolution.y(), false);
  }
#pragma omp parallel for
  // While setting remaining cells to infinity
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++) {
      double cellPhi;
      int cellSign;
      Material::FluidMaterial cellMaterial;
      cellMaterial = _material[i][j];
      cellPhi = _phi[i][j];
      cellSign = cellPhi / std::fabs(cellPhi);

      // Look for surface cells and compute distance over them
      if (i < _resolution.x() - 1 &&
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
      if (j < _resolution.y() - 1 &&
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
    int i = cellId % _resolution.x(), j = cellId / _resolution.x();
    // std::cerr << "Id " << cellId << " -> " << i << ' ' << j << std::endl;
    // If the cell is not processed yet, compute its distance
    if (!processed[i][j]) {
      processed[i][j] = true;
      double distances[2] = {1e8, 1e8};
      int distanceSignal = _phi[i][j] / std::fabs(_phi[i][j]);
      distances[0] =
          std::min(std::fabs(tempPhi[std::max(0, i - 1)][j]),
                   std::fabs(tempPhi[std::min(_resolution.x() - 1, i + 1)][j]));
      distances[1] =
          std::min(std::fabs(tempPhi[i][std::max(0, j - 1)]),
                   std::fabs(tempPhi[i][std::min(_resolution.y() - 1, j + 1)]));
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
    if (i < _resolution.x() - 1 &&
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
    if (j < _resolution.y() - 1 &&
        cellsAdded.find(ijToId(i, j + 1)) == cellsAdded.end()) {
      cellsQueue.push(
          std::make_pair(std::fabs(_phi[i][j + 1]), ijToId(i, j + 1)));
      cellsAdded.insert(ijToId(i, j + 1));
    }
  }
  // Solve Eikonal function
  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      _phi[i][j] = tempPhi[i][j];
    }
  }
}

void LevelSet2::printLevelSetValue() {
  for (int j = _resolution.y() - 1; j >= 0; j--) {
    for (int i = 0; i < _resolution.x(); i++) {
      std::cout << _phi[i][j] << ' ';
    }
    std::cout << std::endl;
  }
}

std::vector<double> &LevelSet2::operator[](const int i) { return _phi[i]; }

void LevelSet2::solvePressure() {
  // Compute velocity divergent over cell center
  // std::vector<double> divergent(cellCount(), 0.0);
  int nCells = cellCount();
  Eigen::VectorXd divergent, pressure;
  Eigen::SparseMatrix<double> pressureMatrix(nCells + 1, nCells + 1);
  std::vector<Eigen::Triplet<double>> triplets;

  divergent = Eigen::VectorXd::Zero(nCells + 1);
  pressure = Eigen::VectorXd::Zero(nCells + 1);

// Solve pressure Poisson equation
#pragma omp parallel for
  // Assemble Laplacian matrix
  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      if (_material[i][j] != Material::FluidMaterial::FLUID)
        continue;

      int validCells = 0;
      int id = ijToId(i, j);
      std::vector<Eigen::Triplet<double>> threadTriplet;

      // TODO: impose second order boundary condition for pressure
      if (i > 0) {
        validCells++;
        if (_material[i - 1][j] != Material::FluidMaterial::AIR)
          threadTriplet.emplace_back(id, ijToId(i - 1, j),
                                     -1 / (_h.x() * _h.x()));
        else {
          double theta = _phi[i][j] / (_phi[i][j] - _phi[i - 1][j]);
        }
      }
      if (i < _resolution.x() - 1) {
        validCells++;
        if (_material[i + 1][j] != Material::FluidMaterial::AIR)
          threadTriplet.emplace_back(id, ijToId(i + 1, j),
                                     -1 / (_h.x() * _h.x()));
        else {
          double theta = _phi[i][j] / (_phi[i][j] - _phi[i + 1][j]);
        }
      }
      if (j > 0) {
        validCells++;
        if (_material[i][j - 1] != Material::FluidMaterial::AIR)
          threadTriplet.emplace_back(id, ijToId(i, j - 1),
                                     -1 / (_h.x() * _h.x()));
        else {
          double theta = _phi[i][j] / (_phi[i][j] - _phi[i][j - 1]);
        }
      }
      if (j < _resolution.y() - 1) {
        validCells++;
        if (_material[i][j + 1] != Material::FluidMaterial::AIR)
          threadTriplet.emplace_back(id, ijToId(i, j + 1),
                                     -1 / (_h.x() * _h.x()));
        else {
          double theta = _phi[i][j] / (_phi[i][j] - _phi[i][j + 1]);
        }
      }
      threadTriplet.emplace_back(id, ijToId(i, j),
                                 validCells / (_h.x() * _h.x()));
#pragma omp critical
      {
        triplets.insert(triplets.end(), threadTriplet.begin(),
                        threadTriplet.end());
      }

      divergent[ijToId(i, j)] = 0;
      divergent[ijToId(i, j)] -= (_u[i + 1][j].x() - _u[i][j].x()) / _h.x();
      divergent[ijToId(i, j)] -= (_v[i][j + 1].y() - _v[i][j].y()) / _h.y();
      divergent[id] /= _dt;
    }
  }
  triplets.emplace_back(nCells, nCells, 1);
  pressureMatrix.setFromTriplets(triplets.begin(), triplets.end());

  // SOlve pressure Poisson system
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;
  solver.compute(pressureMatrix);
  pressure = solver.solve(divergent);
  // std::cerr << divergent.transpose() << std::endl;
  // std::cerr << pressure.transpose() << std::endl;
  // Correct velocity through pressure gradient
  for (int j = 0; j < _resolution.y(); j++) {
    for (int i = 1; i < _resolution.x(); i++) {
      _u[i][j].x(_u[i][j].x() -
                 _dt * (pressure[ijToId(i, j)] - pressure[ijToId(i - 1, j)]) /
                     _h.x());
    }
  }

  for (int j = 1; j < _resolution.y(); j++) {
    for (int i = 0; i < _resolution.x(); i++) {
      _v[i][j].y(_v[i][j].y() -
                 _dt * (pressure[ijToId(i, j)] - pressure[ijToId(i, j - 1)]) /
                     _h.y());
    }
  }
}

} // namespace Ramuh