#include <geometry/bounding_box.h>

namespace Ramuh {

BoundingBox3::BoundingBox3() {}

BoundingBox3::BoundingBox3(Eigen::Array3d min, Eigen::Array3d max) {
  _min = min;
  _max = max;
}

Eigen::Array3d BoundingBox3::min() { return _min; }

Eigen::Array3d BoundingBox3::max() { return _max; }

Eigen::Array3d BoundingBox3::size() { return _max - _min; }

BoundingBox3 BoundingBox3::unitBox() {
  return BoundingBox3(Eigen::Array3d(0, 0, 0), Eigen::Array3d(1, 1, 1));
}

} // namespace Ramuh