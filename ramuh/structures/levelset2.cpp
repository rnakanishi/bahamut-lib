#include <structures/levelset2.h>
#include <utils/macros.h>
#include <omp.h>
#include <cmath>

namespace Ramuh {

LevelSet2::LevelSet2() : RegularGrid2() {
  _phi.resize(_resolution.x());
  for (auto &row : _phi) {
    row.resize(_resolution.y());
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
    row.resize(_resolution.y());
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
      _phi[i][j] =
          std::min(_phi[i][j], position.dot(position) - radius * radius);
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
      if (_phi[i][j] == 0) {
        std::cerr << '[' << i << ',' << j << ']';
        std::cerr << position << distanceLower << distanceUpper << std::endl;
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

std::vector<double> &LevelSet2::operator[](const int i) { return _phi[i]; }

} // namespace Ramuh