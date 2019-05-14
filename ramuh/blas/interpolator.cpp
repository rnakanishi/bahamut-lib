#include <blas/interpolator.h>

namespace Ramuh {

Interpolator::Interpolator() {}

double Interpolator::linear(double target, std::vector<double> points,
                            std::vector<double> values) {
  // TODO: This function assumes that only two points are passed as parameter
  double position = target - points[0];
  position /= (points[1] - points[0]);
  double finalValue = position * (values[1] - values[0]);
  finalValue += values[0];
  return finalValue;
}

} // namespace Ramuh