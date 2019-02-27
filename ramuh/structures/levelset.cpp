#include <structures/levelset.h>
#include <utils/macros.h>
#include <omp.h>

namespace Ramuh {

LevelSet::LevelSet() : RegularGrid() {
  _phi.resize(_resolution.x() + 1);
  for (auto &row : _phi) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1, 1e8);
  }
  _gradPhi.resize(_resolution.x() + 1);
  for (auto &row : _gradPhi) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
  }
  _dt = 0.05;
}

void LevelSet::setResolution(Vector3i newResolution) {
  RegularGrid::setResolution(newResolution);
  _phi.resize(_resolution.x() + 1);
  for (auto &row : _phi) {
    row.resize(_resolution.y() + 1);
    for (auto &depth : row)
      depth.resize(_resolution.z() + 1);
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

void LevelSet::printVertexVelocity() {
  std::cerr << "==== x component: \n";
  for (int k = 0; k < _resolution.z() + 1; k++) {
    for (int j = 0; j < _resolution.y() + 1; j++) {
      for (int i = 0; i < _resolution.x() + 1; i++) {
        std::cerr << _velocity[i][j][k].x() << " ";
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }
  std::cerr << "==== y component: \n";
  for (int k = 0; k < _resolution.z() + 1; k++) {
    for (int j = 0; j < _resolution.y() + 1; j++) {
      for (int i = 0; i < _resolution.x() + 1; i++) {
        std::cerr << _velocity[i][j][k].y() << " ";
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }
}

void LevelSet::addSphereSurface(Vector3d center, double radius) {
  for (int k = 0; k < _resolution.z() + 1; k++)
    for (int j = 0; j < _resolution.y() + 1; j++)
      for (int i = 0; i < _resolution.x() + 1; i++) {
        Vector3d position(Vector3i(i, j, k) * _h);
        position = position - center;
        _phi[i][j][k] =
            std::min(_phi[i][j][k], position.dot(position) - radius * radius);
      }
  checkCellMaterial();
}

void LevelSet::checkCellMaterial() {
  for (int k = 0; k < _resolution.z(); k++)
    for (int j = 0; j < _resolution.y(); j++)
      for (int i = 0; i < _resolution.x(); i++) {
        if (_phi[i][j][k] <= 0 || _phi[i + 1][j][k] <= 0 ||
            _phi[i][j + 1][k] <= 0 || _phi[i][j][k + 1] <= 0 ||
            _phi[i + 1][j + 1][k] <= 0 || _phi[i + 1][j][k + 1] <= 0 ||
            _phi[i][j + 1][k + 1] <= 0 || _phi[i + 1][j + 1][k + 1] <= 0)
          _material[i][j][k] = Material::FluidMaterial::FLUID;
        else
          _material[i][j][k] = Material::FluidMaterial::AIR;
      }
}

void LevelSet::addImplicitFunction() { NOT_IMPLEMENTED(); }

void LevelSet::interpolateVelocitiesToVertices() {

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

void LevelSet::integrateLevelSet() {
  int cellCount;
  cellCount = _resolution.x() * _resolution.y() * _resolution.z();

#pragma omp parallel
  {
    // Compute levelset gradient values
    // Upwind
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    for (int i = id; i < _resolution.x() + 1; i += nthreads)
      for (int j = 0; j < _resolution.y() + 1; j++)
        for (int k = 0; k < _resolution.z() + 1; k++) {
          double dx, dy, dz;
          // X Component
          if (i < _resolution.x())
            dx = (_phi[i + 1][j][k] - _phi[i][j][k]) / _h.x();
          else
            dx = (_phi[i][j][k] - _phi[i - 1][j][k]) / _h.x();

          // Y Component
          if (j < _resolution.y())
            dy = (_phi[i][j + 1][k] - _phi[i][j][k]) / _h.y();
          else
            dy = (_phi[i][j][k] - _phi[i][j - 1][k]) / _h.y();

          // Z Component
          if (k < _resolution.z())
            dz = (_phi[i][j][k + 1] - _phi[i][j][k]) / _h.z();
          else
            dz = (_phi[i][j][k] - _phi[i][j][k - 1]) / _h.z();

          _gradPhi[i][j][k] = Vector3d(dx, dy, dz);
        }
  }

#pragma omp parallel
  {
    // Euler method integration
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    for (int i = id; i < _resolution.x() + 1; i += nthreads)
      for (int j = 0; j < _resolution.y() + 1; j++)
        for (int k = 0; k < _resolution.z() + 1; k++) {
          _phi[i][j][k] -= (_dt * _gradPhi[i][j][k].dot(_velocity[i][j][k]));
        }
  }
}

std::vector<std::vector<double>> &LevelSet::operator[](const int i) {
  return _phi[i];
}

} // namespace Ramuh