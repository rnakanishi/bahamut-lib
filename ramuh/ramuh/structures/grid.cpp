#include <structures/grid.h>

namespace Ramuh {

RegularGrid::RegularGrid() {
  _resolution = Vector3i(32, 32, 32);
  _h = Vector3d(1. / 32);
  _domainSize = Vector3d(1.0);
}

Vector3d RegularGrid::domainSize() { return _domainSize; }

Vector3i RegularGrid::resolution() { return _resolution; }

Vector3d RegularGrid::h() { return _h; }

Vector3d RegularGrid::gridSize() { return _domainSize; }

void RegularGrid::setSize(Vector3d newSize) {
  _domainSize = newSize;
  setH(_domainSize / _resolution);
}

void RegularGrid::setResolution(Vector3i newResolution) {
  _resolution = newResolution;
  setH(_domainSize / _resolution);
  _u.resize(_resolution.x());
  _v.resize(_resolution.y());
  _w.resize(_resolution.z());
}

void RegularGrid::setH(Vector3d newH) { _h = newH; }

} // namespace Ramuh