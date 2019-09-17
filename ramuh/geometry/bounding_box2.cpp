#include <geometry/bounding_box.h>

namespace Ramuh {

BoundingBox2::BoundingBox2() {
  BoundingBox2(Eigen::Array2d(-1), Eigen::Array2d(1));
}

BoundingBox2::BoundingBox2(double min, double max)
    : BoundingBox2(Eigen::Array2d(min), Eigen::Array2d(max)) {}

BoundingBox2::BoundingBox2(Eigen::Array2d min, Eigen::Array2d max) {
  _min = min;
  _max = max;
}

Eigen::Array2d BoundingBox2::min() { return _min; }

Eigen::Array2d BoundingBox2::max() { return _max; }

Eigen::Array2d BoundingBox2::center() { return (_max + _min) / 2.0; }

Eigen::Array2d BoundingBox2::size() { return _max - _min; }

Eigen::Array2d BoundingBox2::clamp(Eigen::Array2d point) {
  point[0] = std::min(_max[0], std::max(_min[0], point[0]));
  point[1] = std::min(_max[1], std::max(_min[1], point[1]));
  return point;
}

bool BoundingBox2::contains(Eigen::Array2d point) {
  if ((point[0] <= _max[0] && point[1] <= _max[1]) &&
      (point[0] >= _min[0] && point[1] >= _min[1]))
    return true;
  return false;
}

bool BoundingBox2::contains(BoundingBox2 box) {
  auto bmin = box.min();
  auto bmax = box.max();
  if ((bmin[0] < _min[0] || bmin[1] < _min[1]) ||
      (bmax[0] > _max[0] || bmax[1] > _max[1]))
    return false;
  return true;
}

BoundingBox2 BoundingBox2::unitBox() {
  return BoundingBox2(Eigen::Array2d(0, 0), Eigen::Array2d(1, 1));
}

} // namespace Ramuh