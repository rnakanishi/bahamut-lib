#include <blas/interpolator.h>

namespace Ramuh {

Interpolator::Interpolator() {}

double Interpolator::linear(double target, std::vector<double> samples,
                            std::vector<double> values) {
  // TODO: This function assumes that only two samples are passed as parameter
  double position = target - samples[0];
  position /= (samples[1] - samples[0]);
  double finalValue = position * (values[1] - values[0]);
  finalValue += values[0];
  return finalValue;
}

double Interpolator::trilinear(Eigen::Array3d position,
                               std::vector<Eigen::Array3d> samples,
                               std::vector<double> sampleValues) {
  double intermediateValues[4];
  for (int it = 0; it < 4; it++) {
    std::vector<double> points, values;
    points.push_back(samples[2 * it + 0][0]);
    points.push_back(samples[2 * it + 1][0]);

    values.push_back(sampleValues[2 * it + 0]);
    values.push_back(sampleValues[2 * it + 1]);

    intermediateValues[it] = Interpolator::linear(position[0], points, values);
  }

  for (int it = 0; it < 2; it++) {
    std::vector<double> points, values;
    points.push_back(samples[4 * it + 0][1]);
    points.push_back(samples[4 * it + 2][1]);

    values.push_back(intermediateValues[4 * it + 0]);
    values.push_back(intermediateValues[4 * it + 2]);
    intermediateValues[it] = Interpolator::linear(position[1], points, values);
  }

  double interpolated;
  {
    int it = 0;
    std::vector<double> points, values;
    points.push_back(samples[8 * it + 0][2]);
    points.push_back(samples[8 * it + 4][2]);

    values.push_back(intermediateValues[0]);
    values.push_back(intermediateValues[1]);
    interpolated = Interpolator::linear(position[2], points, values);
  }
  return interpolated;
}
} // namespace Ramuh