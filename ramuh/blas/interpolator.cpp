#include <blas/interpolator.h>
#include <iostream>
#include <Eigen/Dense>

namespace Ramuh {

Interpolator::Interpolator() {}

double Interpolator::linear(double target, std::vector<double> points,
                            std::vector<double> values) {
  // TODO: This function assumes that only two points are passed as parameter
  double theta =
      (target - points[0]) * (values[1] - values[0]) / (points[1] - points[0]);
  return theta + values[0];
}

double Interpolator::bilinear(double position[2],
                              std::vector<Eigen::Array2d> samples,
                              std::vector<double> sampleValues) {
  double intermediateValues[2];
  for (int it = 0; it < 2; it++) {
    std::vector<double> points, values;
    points.push_back(samples[2 * it + 0][0]);
    points.push_back(samples[2 * it + 1][0]);

    values.push_back(sampleValues[2 * it + 0]);
    values.push_back(sampleValues[2 * it + 1]);
    intermediateValues[it] = Interpolator::linear(position[0], points, values);
  }

  double interpolated;
  for (int it = 0; it < 1; it++) {
    std::vector<double> points, values;
    points.push_back(samples[4 * it + 0][1]);
    points.push_back(samples[4 * it + 2][1]);

    values.push_back(intermediateValues[2 * it + 0]);
    values.push_back(intermediateValues[2 * it + 1]);
    interpolated = Interpolator::linear(position[1], points, values);
  }
  return interpolated;
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

    values.push_back(intermediateValues[2 * it + 0]);
    values.push_back(intermediateValues[2 * it + 1]);
    intermediateValues[it] = Interpolator::linear(position[1], points, values);
  }

  double interpolated;
  {
    int it = 0;
    std::vector<double> points, values;
    points.push_back(samples[8 * it + 0][2]);
    points.push_back(samples[8 * it + 4][2]);

    values.push_back(intermediateValues[2 * it + 0]);
    values.push_back(intermediateValues[2 * it + 1]);
    interpolated = Interpolator::linear(position[2], points, values);
  }
  return interpolated;
}
} // namespace Ramuh