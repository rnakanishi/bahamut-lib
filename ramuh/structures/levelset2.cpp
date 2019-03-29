#include <structures/levelset2.h>
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
      position = position - center;
      double value = position.dot(position) - radius * radius;
      int sign = value / std::fabs(value);
      _phi[i][j] = sign * std::min(std::fabs(_phi[i][j]), std::fabs(value));
    }
  checkCellMaterial();
}

void LevelSet2::addCubeSurface(Vector2d lower, Vector2d upper) {
  // for (int k = 0; k < _resolution.z(); k++)
  for (int j = 0; j < _resolution.y(); j++)
    for (int i = 0; i < _resolution.x(); i++) {
      Vector2d position(Vector2i(i, j) * _h + _h / 2);
      Vector2d distanceLower = (position - lower).abs();
      Vector2d distanceUpper = (position - upper).abs();
      if (position > lower && position < upper) {
        _phi[i][j] = std::min(
            _phi[i][j], -std::min(distanceLower.min(), distanceUpper.min()));
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

#pragma omp for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++) {
      Vector2d position, backPosition, velocity, h(_h.x(), _h.y());
      Vector2i index;
      double newPhi = 0.0;
      double distanceCount = 0.0, distance = 0.0;

      position = Vector2d(i, j) * h + h / 2.0;
      velocity.x((_u[i][j].x() + _u[i + 1][j].x()) / 2.0);
      velocity.y((_v[i][j].y() + _v[i][j + 1].y()) / 2.0);
      backPosition = position - velocity * _dt;

      index.x(std::floor(backPosition.x() / h.x()));
      index.y(std::floor(backPosition.y() / h.y()));

      if (index >= Vector2i(0, 0) &&
          index < Vector2(_resolution.x(), _resolution.y())) {
        Material::FluidMaterial centerMaterial =
            _material[index.x()][index.y()];
        position = index * h + h / 2.0;
        distance = (backPosition - position).length();
        distanceCount += distance;
        newPhi += distance * oldPhi[index.x()][index.y()];
        if (index.x() > 0) {
          //  && _material[index.x() - 1][index.y()][0] == centerMaterial) {
          distance = (backPosition - position).length();
          distanceCount += distance;
          newPhi += distance * oldPhi[index.x() - 1][index.y()];
        }
        if (index.x() < _resolution.x() - 1) {
          //  && _material[index.x() + 1][index.y()][0] == centerMaterial) {
          distance = (backPosition - position).length();
          distanceCount += distance;
          newPhi += distance * oldPhi[index.x() + 1][index.y()];
        }
        if (index.y() > 0) {
          // && _material[index.x()][index.y() - 1][0] == centerMaterial) {
          distance = (backPosition - position).length();
          distanceCount += distance;
          newPhi += distance * oldPhi[index.x()][index.y() - 1];
        }
        if (index.y() < _resolution.y() - 1) {
          // && _material[index.x()][index.y() + 1][0] == centerMaterial) {
          distance = (backPosition - position).length();
          distanceCount += distance;
          newPhi += distance * oldPhi[index.x()][index.y() + 1];
        }
        // TODO: Do the same for z axis
      } else {
        newPhi = oldPhi[i][j];
        distanceCount += 1.0;
      }
      if (distanceCount == 0.0) {
        _phi[i][j] = oldPhi[i][j];
      } else
        _phi[i][j] = newPhi / distanceCount;
    }

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
      cellSign = (cellPhi < 0.0) ? -1 : 1;

      // Look for surface cells and compute distance over them
      if (i < _resolution.x() - 1 &&
          std::signbit(cellPhi) != std::signbit(_phi[i + 1][j])) {
        // compute distance to the surface
        double distance = 1e8;
        double theta = cellSign * _phi[i][j] / (_phi[i][j] - _phi[i + 1][j]);
        tempPhi[i][j] = std::min(tempPhi[i][j], theta * _h.x());
        processed[i][j] = true;
      }
      if (i > 0 && std::signbit(cellPhi) != std::signbit(_phi[i - 1][j])) {
        // compute distance to the surface
        double distance = 1e8;
        double theta = cellSign * _phi[i][j] / (_phi[i][j] - _phi[i - 1][j]);
        tempPhi[i][j] = std::min(tempPhi[i][j], theta * _h.x());
        processed[i][j] = true;
      }
      if (j < _resolution.y() - 1 &&
          std::signbit(cellPhi) != std::signbit(_phi[i][j + 1])) {
        // compute distance to the surface
        double distance = 1e8;
        double theta = cellSign * _phi[i][j] / (_phi[i][j] - _phi[i][j + 1]);
        tempPhi[i][j] = std::min(tempPhi[i][j], theta * _h.x());
        processed[i][j] = true;
      }
      if (j > 0 && std::signbit(cellPhi) != std::signbit(_phi[i][j - 1])) {
        // compute distance to the surface
        double distance = 1e8;
        double theta = cellSign * _phi[i][j] / (_phi[i][j] - _phi[i][j - 1]);
        tempPhi[i][j] = std::min(tempPhi[i][j], theta * _h.x());
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
      distances[0] = std::min(tempPhi[std::max(0, i - 1)][j],
                              tempPhi[std::min(_resolution.x() - 1, i + 1)][j]);
      distances[1] = std::min(tempPhi[i][std::max(0, j - 1)],
                              tempPhi[i][std::min(_resolution.y() - 1, j + 1)]);
      if (std::fabs(distances[0]) > std::fabs(distances[1])) {
        double aux = distances[0];
        distances[0] = distances[1];
        distances[1] = aux;
      }
      // std::cerr << "Distances " << distances[0] << ' ' << distances[1]
      // << std::endl;
      double newPhi = distances[0] + ((distances[0] < 0) ? -_h.x() : _h.x());
      if (newPhi > distances[1]) {
        newPhi =
            0.5 *
            (distances[0] + distances[1] +
             std::sqrt(2 * _h.x() * _h.x() - ((distances[1] - distances[0]) *
                                              (distances[1] - distances[0]))));
      }
      if (newPhi < tempPhi[i][j])
        tempPhi[i][j] = newPhi;
    }
    // Compute neighbor distance to the levelset and add it to queue
    std::cerr << cellsAdded.size() << ' ';
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
      std::cerr << tempPhi[i][j] << ' ';
    }
    std::cerr << '\n';
  }
} // namespace Ramuh

std::vector<double> &LevelSet2::operator[](const int i) { return _phi[i]; }

} // namespace Ramuh