#include <geometry/bounding_box.h>
#include <utils/log.hpp>

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

Eigen::Array2d BoundingBox2::getMin() { return _min; }

Eigen::Array2d BoundingBox2::getMax() { return _max; }

Eigen::Array2d BoundingBox2::getCenter() { return (_max + _min) / 2.0; }

Eigen::Array2d BoundingBox2::getSize() { return _max - _min; }

void BoundingBox2::setMin(Eigen::Array2d newMin) { _min = newMin; }

void BoundingBox2::setMax(Eigen::Array2d newMax) { _max = newMax; }

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
  auto bmin = box.getMin();
  auto bmax = box.getMax();
  if ((bmin[0] < _min[0] || bmin[1] < _min[1]) ||
      (bmax[0] > _max[0] || bmax[1] > _max[1]))
    return false;
  return true;
}

BoundingBox2 BoundingBox2::unitBox() {
  return BoundingBox2(Eigen::Array2d(0, 0), Eigen::Array2d(1, 1));
}

Eigen::Array2d BoundingBox2::findIntersection(Eigen::Array2d p1,
                                              Eigen::Array2d p2) {
  if (!contains(p1) && !contains(p2)) {
    Ramuh::Log::raiseError(
        "BoundingBox2::findIntersection: segment does not intersect bbox.");
    return Eigen::Array2d(0, 0);
  }
  Eigen::Array2d origin, ending;
  origin = ending = p1;
  if (contains(p1)) {
    origin = p2;
  } else {
    ending = p2;
  }
  double angle = (p2[1] - p1[1]) / (p2[0] - p1[0]);
  Eigen::Vector2d direction;
  direction = (p2 - p1).matrix().normalized();
  Eigen::Array2d newOrigin = origin;

  // for each edge, check intersection. If true, project to that edge
  // bottom/top
  if (ending[1] < _min[1]) {
    ending[1] = _min[1];
    ending[0] = (ending[1] - newOrigin[1]) / angle + origin[0];
  } else if (ending[1] > _max[1]) {
    ending[1] = _max[1];
    ending[0] = (ending[1] - newOrigin[1]) / angle + origin[0];
  }
  // left/right
  if (ending[0] < _min[0]) {
    ending[0] = _min[0];
    ending[1] = angle * (ending[0] - newOrigin[0]) + newOrigin[1];
  } else if (ending[0] > _max[0]) {
    ending[0] = _max[0];
    ending[1] = angle * (ending[0] - newOrigin[0]) + newOrigin[1];
  }

  return ending;
}

} // namespace Ramuh