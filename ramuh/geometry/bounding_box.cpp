#include <geometry/bounding_box.h>

namespace Ramuh {

BoundingBox1::BoundingBox1() {}

BoundingBox1::BoundingBox1(double min, double max) {
  _min = min;
  _max = max;
}

double BoundingBox1::min() { return _min; }

double BoundingBox1::max() { return _max; }

double BoundingBox1::center() { return (_max + _min) / 2.0; }

double BoundingBox1::size() { return _max - _min; }

double BoundingBox1::clamp(double point) {
  point = std::min(_max, std::max(_min, point));
  return point;
}

bool BoundingBox1::contains(double point) {
  if (point <= _max && point >= _min)
    return true;
  return false;
}

BoundingBox1 BoundingBox1::unitBox() { return BoundingBox1(0., 1.); }

} // namespace Ramuh