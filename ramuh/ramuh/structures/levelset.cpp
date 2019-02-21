#include <structures/levelset.h>
#include <utils/macros.h>
#include <omp.h>

namespace Ramuh {

LevelSet::LevelSet() : RegularGrid() {
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
  dt = 0.05;
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
}

void LevelSet::addImplicitFunction() { NOT_IMPLEMENTED(); }

void LevelSet::interpolateVelocitiesToVertices() { NOT_IMPLEMENTED(); }

void LevelSet::integrateLevelSet() {
  int cellCount;
  cellCount = _resolution.x() * _resolution.y() * _resolution.z();

// Compute levelset gradient values
// Upwind
#pragma omp parallel
  {
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();

    for (int i = id; i < _resolution.x(); i += nthreads) {
      for (int j = 0; j < _resolution.y(); j++) {
        for (int k = 0; k < _resolution.z(); k++) {
          double dx, dy, dz;
          // X Component
          if (i < _resolution.x())
            dx = (_phi[i + 1][j][k] - _phi[i][j][k]) / _h.x();
          else
            dx = (_phi[i][j][k] - _phi[i - 1][j][k]) / _h.x();

          // Y Component
          if (i < _resolution.y())
            dy = (_phi[i][j + 1][k] - _phi[i][j][k]) / _h.y();
          else
            dy = (_phi[i][j][k] - _phi[i][j - 1][k]) / _h.y();

          // Z Component
          if (i < _resolution.z())
            dz = (_phi[i][j][k + 1] - _phi[i][j][k]) / _h.z();
          else
            dz = (_phi[i][j][k] - _phi[i][j][k - 1]) / _h.z();

          _gradPhi[i][j][k] = Vector3d(dx, dy, dz);
        }
      }
    }
  }
  // Euler method integration
}

void LevelSet::solvePressure() {

  // Compute velocity divergent over cell center

  // Solve pressure Poisson equation

  // Correct velocity through pressure gradient
}

std::vector<std::vector<double>> &LevelSet::operator[](const int i) {
  return _phi[i];
}
} // namespace Ramuh