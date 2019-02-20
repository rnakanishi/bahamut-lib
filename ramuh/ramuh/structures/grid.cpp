#include <structures/grid.h>

namespace Ramuh {

RegularGrid::RegularGrid() {
  _resolution = Vector3i(32, 32, 32);
  _h = Vector3d(1. / 32);
  _size = Vector3d(1.0);
}

Vector3d RegularGrid::size() { return _size; }

Vector3i RegularGrid::resolution() { return _resolution; }

Vector3d RegularGrid::h() { return _h; }

Vector3d RegularGrid::gridSize() { return _size; }

void RegularGrid::setSize(Vector3d newSize) {
  _size = newSize;
  setH(_size / _resolution);
}

void RegularGrid::setResolution(Vector3i newResolution) {
  _resolution = newResolution;
  setH(_size / _resolution);
}

void RegularGrid::setH(Vector3d newH) { _h = newH; }

} // namespace Ramuh