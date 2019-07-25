#include <geometry/bounding_box.h>

namespace Ramuh {

BoundingBox3::BoundingBox3() {}

BoundingBox3::BoundingBox3(Eigen::Array3d min, Eigen::Array3d max) {
  _min = min;
  _max = max;
}

Eigen::Array3d BoundingBox3::min() { return _min; }

Eigen::Array3d BoundingBox3::max() { return _max; }

Eigen::Array3d BoundingBox3::center() { return (_max + _min) / 2.0; }

Eigen::Array3d BoundingBox3::size() { return _max - _min; }

Eigen::Array3d BoundingBox3::clamp(Eigen::Array3d point) {
  point[0] = std::min(_max[0], std::max(_min[0], point[0]));
  point[1] = std::min(_max[1], std::max(_min[1], point[1]));
  point[2] = std::min(_max[2], std::max(_min[2], point[2]));
  return point;
}

bool BoundingBox3::contains(Eigen::Array3d point) {
  if ((point[0] <= _max[0] && point[1] <= _max[1] && point[2] <= _max[2]) &&
      (point[0] >= _min[0] && point[1] >= _min[1] && point[2] >= _min[2]))
    return true;
  return false;
}

BoundingBox3 BoundingBox3::unitBox() {
  return BoundingBox3(Eigen::Array3d(0, 0, 0), Eigen::Array3d(1, 1, 1));
}

} // namespace Ramuh