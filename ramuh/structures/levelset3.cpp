#include <structures/levelset3.h>
#include <utils/macros.h>
#include <omp.h>
#include <cmath>

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
  double **oldPhi;
  oldPhi = new double *[_resolution.x()];
  for (int i = 0; i < _resolution.x(); i++) {
    oldPhi[i] = new double[_resolution.y()];
  }

#pragma omp parallel for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++)
      oldPhi[i][j] = _phi[i][j][0];

#pragma omp for
  for (int i = 0; i < _resolution.x(); i++)
    for (int j = 0; j < _resolution.y(); j++) {
      Vector2d position, backPosition, velocity, h(_h.x(), _h.y());
      Vector2i index;
      double newPhi = 0.0;
      double distanceCount = 0.0, distance = 0.0;

      position = Vector2d(i, j) * h + h / 2.0;
      velocity.x((_u[i][j][0].x() + _u[i + 1][j][0].x()) / 2.0);
      velocity.y((_v[i][j][0].y() + _v[i][j + 1][0].y()) / 2.0);
      backPosition = position - velocity * _dt;

      index.x(std::floor(backPosition.x() / h.x()));
      index.y(std::floor(backPosition.y() / h.y()));

      if (index >= Vector2i(0, 0) &&
          index < Vector2(_resolution.x(), _resolution.y())) {
        Material::FluidMaterial centerMaterial =
            _material[index.x()][index.y()][0];
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
        _phi[i][j][0] = oldPhi[i][j];
      } else
        _phi[i][j][0] = newPhi / distanceCount;
    }

  for (int i = 0; i < _resolution.y(); ++i) {
    delete[] oldPhi[i];
  }
  delete[] oldPhi;
}

std::vector<std::vector<double>> &LevelSet3::operator[](const int i) {
  return _phi[i];
}

} // namespace Ramuh