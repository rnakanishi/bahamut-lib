#include <structures/levelset3.h>
#include <geometry/matrix.h>
#include <utils/macros.h>
#include <omp.h>
#include <queue>
#include <cmath>
#include <set>

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
        Vector3d position(Vector3i(i, j, 0) * _h + _h / 2);
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
  int cellCount;
  cellCount = _resolution.x() * _resolution.y();
  // std::vector<std::vector<double>> oldPhi;
  Matrix3<double> oldPhi(_resolution);

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int k = 0; j < _resolution.z(); k++)
        oldPhi[i][j][k] = _phi[i][j][k];

#pragma omp for
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
            iCandidates.push_back(index.x() - 1);
          else if (backPosition.x() < cellCenter.x() && index.x() > 0)
            iCandidates.push_back(index.x() + 1);
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

  tempPhi.changeSize(_resolution);
}

} // namespace Ramuh