#include <structures/levelset3.h>
#include <geometry/matrix.h>
#include <utils/macros.h>
#include <omp.h>
#include <queue>
#include <cmath>
#include <set>
#include <utility>
#include <algorithm>

namespace Ramuh {

LevelSet3::LevelSet3() : RegularGrid3() {
  _phi.resize(_resolution.x());
  for (auto &row : _phi) {
    row.resize(_resolution.y());
    for (auto &depth : row)
      depth.resize(_resolution.z(), 1e6);
  }
  _gradPhi.resize(_resolution.x() + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
  }
}

void LevelSet3::setResolution(Vector3i newResolution) {
  RegularGrid3::setResolution(newResolution);
  _phi.resize(_resolution.x());
  for (auto &row : _phi) {
    row.resize(_resolution.y());
    for (auto &depth : row)
      depth.resize(_resolution.z(), 1e6);
  }
  _gradPhi.resize(_resolution.x() + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
  }
  _velocity.resize(_resolution.x() + 1);
  for (auto &row : _velocity) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
  }
}

void LevelSet3::printVertexVelocity() {
  std::cerr << "==== x component: \n";
  for (int j = 0; j < _resolution.y() + 1; j++) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _velocity[i][j][0].x() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "==== y component: \n";
  for (int j = 0; j < _resolution.y() + 1; j++) {
    for (int i = 0; i < _resolution.x() + 1; i++) {
      std::cerr << _velocity[i][j][0].y() << " ";
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
}

void LevelSet3::addSphereSurface(Vector3d center, double radius) {
  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int i = 0; i < _resolution.x(); i++) {
        Vector3d position(Vector3i(i, j, k) * _h + _h / 2);
        double distance = (position - center).length() - radius;
        _phi[i][j][k] = distance;
      }
  checkCellMaterial();
}

void LevelSet3::addCubeSurface(Vector3d lower, Vector3d upper) {
  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int i = 0; i < _resolution.x(); i++) {
        Vector3d position(Vector3i(i, j, k) * _h + _h / 2);
        Vector3d distanceLower = (position - lower).abs();
        Vector3d distanceUpper = (position - upper).abs();
        if (position > lower && position < upper) {
          _phi[i][j][k] =
              -std::min(_phi[i][j][k],
                        std::min(distanceLower.min(), distanceUpper.min()));
        } else {
          _phi[i][j][k] =
              std::min(_phi[i][j][k],
                       std::min(distanceLower.min(), distanceUpper.min()));
        }
      }
  checkCellMaterial();
}

void LevelSet3::checkCellMaterial() {
  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int i = 0; i < _resolution.x(); i++) {
        if (_phi[i][j][k] <= 0)
          _material[i][j][k] = Material::FluidMaterial::FLUID;
        else
          _material[i][j][k] = Material::FluidMaterial::AIR;
      }
}

void LevelSet3::addImplicitFunction() { NOT_IMPLEMENTED(); }

void LevelSet3::interpolateVelocitiesToVertices() {

// TODO: Improve this interpolation
#pragma omp parallel
  {
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    // Iterate over all grid vertices interpolating faces velocities
    for (int i = id; i < _resolution.x() + 1; i += nthreads)
      for (int j = 0; j < _resolution.y() + 1; j++)
        for (int k = 0; k < _resolution.z() + 1; k++) {
          Vector3d vel;
          // u component
          if (i == 0 || i == _resolution.x())
            vel.x(0.0);
          else {
            if (j > 0 && k > 0 && j < _resolution.y() && k < _resolution.z())
              vel.x((_u[i][j][k].x() + _u[i][j - 1][k].x() +
                     _u[i][j - 1][k - 1].x() + _u[i][j][k - 1].x()) /
                    4.0);
            else if (j == 0 && k == 0)
              vel.x(_u[i][j][k].x());
            else if (j == 0 && k > 0 && k < _resolution.z())
              vel.x((_u[i][j][k].x() + _u[i][j][k - 1].x()) / 2.0);
            else if (j == 0 && k == _resolution.z())
              vel.x(_u[i][j][k - 1].x());
            else if (k == 0 && j > 0 && j < _resolution.y())
              vel.x((_u[i][j][k].x() + _u[i][j - 1][k].x()) / 2.0);
            else if (k == 0 && j == _resolution.y())
              vel.x(_u[i][j - 1][k].x());
            else if (j == _resolution.y() && k > 0 && k < _resolution.z())
              vel.x((_u[i][j - 1][k].x() + _u[i][j - 1][k - 1].x()) / 2.0);
            else if (k == _resolution.z() && j > 0 && j < _resolution.y())
              vel.x((_u[i][j - 1][k - 1].x() + _u[i][j][k - 1].x()) / 2.0);
            else
              vel.x(_u[i][j - 1][k - 1].x());
          }
          // v component
          if (j == 0 || j == _resolution.y())
            vel.y(0.0);
          else {
            if (i > 0 && k > 0 && i < _resolution.x() && k < _resolution.z())
              vel.y((_v[i][j][k].y() + _v[i - 1][j][k].y() +
                     _v[i - 1][j][k - 1].y() + _v[i][j][k - 1].y()) /
                    4.0);
            else if (i == 0 && k == 0)
              vel.y(_v[i][j][k].y());
            else if (i == 0 && k > 0 && k < _resolution.z())
              vel.y((_v[i][j][k].y() + _v[i][j][k - 1].y()) / 2.0);
            else if (i == 0 && k == _resolution.z())
              vel.y((_v[i][j][k - 1].y()));
            else if (k == 0 && i > 0 && i < _resolution.x())
              vel.y((_v[i][j][k].y() + _v[i - 1][j][k].y()) / 2.0);
            else if (k == 0 && i == _resolution.x())
              vel.y((_v[i - 1][j][k].y()));
            else if (i == _resolution.x() && k > 0 && k < _resolution.z())
              vel.y((_v[i - 1][j][k].y() + _v[i - 1][j][k - 1].y()) / 2.0);
            else if (k == _resolution.z() && i > 0 && i < _resolution.x())
              vel.y((_v[i - 1][j][k - 1].y() + _v[i][j][k - 1].y()) / 2.0);
            else
              vel.y(_v[i - 1][j][k - 1].y());
          }
          // z component
          if (k == 0 || k == _resolution.z())
            vel.z(0);
          else {
            if (i > 0 && j > 0 && i < _resolution.x() && j < _resolution.y())
              vel.z((_w[i][j][k].z() + _w[i - 1][j][k].z() +
                     _w[i - 1][j - 1][k].z() + _w[i][j - 1][k].z()) /
                    4.0);
            else if (i == 0 && j == 0)
              vel.z(_w[i][j][k].z());
            else if (i == 0 && j > 0 && j < _resolution.y())
              vel.z((_w[i][j][k].z() + _w[i][j - 1][k].z()) / 2.0);
            else if (i == 0 && j == _resolution.y())
              vel.z(_w[i][j - 1][k].z());
            else if (j == 0 && i > 0 && i < _resolution.x())
              vel.z((_w[i][j][k].z() + _w[i - 1][j][k].z()) / 2.0);
            else if (j == 0 && i == _resolution.x())
              vel.z(_w[i - 1][j][k].z());
            else if (i == _resolution.x() && j > 0 && j < _resolution.y())
              vel.z((_w[i - 1][j][k].z() + _w[i - 1][j - 1][k].z()) / 2.0);
            else if (j == _resolution.y() && i > 0 && i < _resolution.x())
              vel.z((_w[i - 1][j - 1][k].z() + _w[i][j - 1][k].z()) / 2.0);
            else
              vel.z(_w[i - 1][j - 1][k].z());
          }

          _velocity[i][j][k] = vel;
        }
  }

} // namespace Ramuh

void LevelSet3::integrateLevelSet() {
  // std::vector<std::vector<double>> oldPhi;
  Matrix3<double> oldPhi;
  oldPhi.changeSize(_resolution);

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++)
        oldPhi[i][j][k] = _phi[i][j][k];

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        Vector3d position, backPosition, velocity, h(_h.x(), _h.y(), _h.z()),
            cellCenter;
        Vector3i index;
        double newPhi = oldPhi[i][j][k];
        double distanceCount = 0.0, distance = 0.0;

        position = Vector3d(i, j, k) * h + h / 2.0;
        velocity.x((_u[i][j][k].x() + _u[i + 1][j][k].x()) / 2.0);
        velocity.y((_v[i][j][k].y() + _v[i][j + 1][k].y()) / 2.0);
        velocity.z((_w[i][j][k].z() + _w[i][j][k + 1].z()) / 2.0);
        backPosition = position - (velocity * _dt);
        index.x(std::floor(backPosition.x() / h.x()));
        index.y(std::floor(backPosition.y() / h.y()));
        index.z(std::floor(backPosition.z() / h.z()));

        if (velocity.length() > 1e-8 && index >= Vector3i(0) &&
            index < _resolution) {
          cellCenter = Vector3d(index.x(), index.y(), index.z()) * h + h / 2.0;
          std::vector<int> iCandidates, jCandidates, kCandidates;
          iCandidates.push_back(index.x());
          jCandidates.push_back(index.y());
          kCandidates.push_back(index.z());
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
          if (backPosition.z() > cellCenter.z() &&
              index.z() < _resolution.z() - 1)
            kCandidates.push_back(index.z() + 1);
          else if (backPosition.z() < cellCenter.z() && index.z() > 0)
            kCandidates.push_back(index.z() - 1);

          newPhi = 0.;
          for (auto u : iCandidates)
            for (auto v : jCandidates)
              for (auto w : kCandidates) {
                position = Vector3d(u, v, w) * h + h / 2.0;
                distance = (backPosition - position).length();
                distanceCount += 1. / distance;
                newPhi += (1. / distance) * oldPhi[u][v][w];
              }
          newPhi /= distanceCount;
        }
        _phi[i][j][k] = newPhi;
      }
}

std::vector<std::vector<double>> &LevelSet3::operator[](const int i) {
  return _phi[i];
}

void LevelSet3::redistance() {
  Matrix3<double> tempPhi;
  Matrix3<bool> processed;
  std::priority_queue<std::pair<double, int>,
                      std::vector<std::pair<double, int>>,
                      std::greater<std::pair<double, int>>>
      cellsQueue;
  std::set<int> cellsAdded;

  tempPhi.changeSize(_resolution, 1e8);
  processed.changeSize(_resolution, false);

#pragma omp parallel for
  // Find surface cells adding them to queue
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++) {
        double cellPhi = _phi[i][j][k];
        int cellSign = cellPhi / std::fabs(cellPhi);

        // For each surface pair, compute their distance and add to queue
        double theta = 1e8;
        if (i < _resolution.x() - 1 &&
            std::signbit(cellPhi) != std::signbit(_phi[i + 1][j][k]))
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i + 1][j][k])));
        if (i > 0 && std::signbit(cellPhi) != std::signbit(_phi[i - 1][j][k]))
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i - 1][j][k])));
        if (j < _resolution.y() - 1 &&
            std::signbit(cellPhi) != std::signbit(_phi[i][j + 1][k]))
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i][j + 1][k])));
        if (j > 0 && std::signbit(cellPhi) != std::signbit(_phi[i][j - 1][k]))
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i][j - 1][k])));
        if (k < _resolution.z() - 1 &&
            std::signbit(cellPhi) != std::signbit(_phi[i][j][k + 1]))
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i][j][k + 1])));
        if (k > 0 && std::signbit(cellPhi) != std::signbit(_phi[i][j][k - 1]))
          theta = std::min(theta,
                           std::fabs(cellPhi / (cellPhi - _phi[i][j][k - 1])));

        if (theta < 1e7) {
          tempPhi[i][j][k] =
              cellSign * std::min(std::fabs(tempPhi[i][j][k]), theta * _h.x());
          processed[i][j][k] = true;
        }

        if (processed[i][j][k]) {
#pragma omp critical
          {
            cellsQueue.push(
                std::make_pair(std::fabs(tempPhi[i][j][k]), ijkToId(i, j, k)));
            cellsAdded.insert(ijkToId(i, j, k));
          }
        }
      }

  // Propagating distances
  while (!cellsQueue.empty()) {
    int cellId = cellsQueue.top().second;
    cellsQueue.pop();
    int k = cellId / (_resolution.x() * _resolution.y());
    int j = (cellId % (_resolution.x() * _resolution.y())) / _resolution.x();
    int i = (cellId % (_resolution.x() * _resolution.y())) % _resolution.x();
    // TODO: check this caculation

    if (!processed[i][j][k]) {
      processed[i][j][k] = true;
      std::vector<double> distances;
      int distanceSignal = _phi[i][j][k] / std::fabs(_phi[i][j][k]);
      distances[0] = std::min(
          std::fabs(tempPhi[std::max(0, i - 1)][j][k]),
          std::fabs(tempPhi[std::min(_resolution.x() - 1, i + 1)][j][k]));
      distances[1] = std::min(
          std::fabs(tempPhi[i][std::max(0, j - 1)][k]),
          std::fabs(tempPhi[i][std::min(_resolution.y() - 1, j + 1)][k]));
      distances[2] = std::min(
          std::fabs(tempPhi[i][j][std::max(0, k - 1)]),
          std::fabs(tempPhi[i][j][std::min(_resolution.z() - 1, k + 1)]));
      std::sort(distances.begin(), distances.end());

      // Checking each component and redistancing
      double newPhi = distances[0] + _h.x();
      if (newPhi > distances[1]) {
        double h = _h.x() * _h.x();
        newPhi =
            0.5 *
            (distances[0] + distances[1] +
             std::sqrt(2 * _h.x() * _h.x() - ((distances[1] - distances[0]) *
                                              (distances[1] - distances[0]))));
        if (std::fabs(newPhi) > distances[2]) {
          // TODO: check calculations from Bridson's book
          newPhi = 0.0;
        }
      }

      if (newPhi < std::fabs(tempPhi[i][j][k]))
        tempPhi[i][j][k] = distanceSignal * newPhi;

      // Add all neighbors of the cell to queue
      if (i > 0 && cellsAdded.find(ijkToId(i - 1, j, k)) == cellsAdded.end()) {
        cellsQueue.push(
            std::make_pair(std::fabs(_phi[i - 1][j][k]), ijkToId(i - 1, j, k)));
        cellsAdded.insert(ijkToId(i - 1, j, k));
      }
      if (i < _resolution.x() - 1 &&
          cellsAdded.find(ijkToId(i + 1, j, k)) == cellsAdded.end()) {
        cellsQueue.push(
            std::make_pair(std::fabs(_phi[i + 1][j][k]), ijkToId(i + 1, j, k)));
        cellsAdded.insert(ijkToId(i + 1, j, k));
      }
      if (j > 0 && cellsAdded.find(ijkToId(i, j - 1, k)) == cellsAdded.end()) {
        cellsQueue.push(
            std::make_pair(std::fabs(_phi[i][j - 1][k]), ijkToId(i, j - 1, k)));
        cellsAdded.insert(ijkToId(i, j - 1, k));
      }
      if (j < _resolution.y() - 1 &&
          cellsAdded.find(ijkToId(i, j + 1, k)) == cellsAdded.end()) {
        cellsQueue.push(
            std::make_pair(std::fabs(_phi[i][j + 1][k]), ijkToId(i, j + 1, k)));
        cellsAdded.insert(ijkToId(i, j + 1, k));
      }
      if (k > 0 && cellsAdded.find(ijkToId(i, j, k - 1)) == cellsAdded.end()) {
        cellsQueue.push(
            std::make_pair(std::fabs(_phi[i][j][k - 1]), ijkToId(i, j, k - 1)));
        cellsAdded.insert(ijkToId(i, j, k - 1));
      }
      if (k < _resolution.z() - 1 &&
          cellsAdded.find(ijkToId(i, j, k + 1)) == cellsAdded.end()) {
        cellsQueue.push(
            std::make_pair(std::fabs(_phi[i][j][k + 1]), ijkToId(i, j, k + 1)));
        cellsAdded.insert(ijkToId(i, j, k + 1));
      }
    }
  }
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; k < _resolution.z(); k++)
        _phi[i][j][k] = tempPhi[i][j][k];
}

void LevelSet3::printLevelSetValue() {

  for (int k = 0; k < _resolution.z(); k++) {
    for (int j = _resolution.y() - 1; j >= 0; j--) {
      for (int i = 0; i < _resolution.x(); i++) {
        std::cout << _phi[i][j][k] << ' ';
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

} // namespace Ramuh