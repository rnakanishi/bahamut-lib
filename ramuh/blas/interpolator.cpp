#include <blas/interpolator.h>
#include <iostream>
#include <Eigen/Dense>
#include <functional>
#include <cmath>

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

double Interpolator::catmullRom(double position, std::vector<double> points,
                                std::vector<double> values) {

  auto tj = [](double t, double points[2], double values[2]) {
    double distance;
    distance = std::sqrt((points[1] - points[0]) * (points[1] - points[0]) +
                         (values[1] - values[0]) * (values[1] - values[0]));
    return t + std::pow(distance, 0.5);
  };

  double t[4];
  t[0] = 0;
  for (size_t i = 1; i < 4; i++) {
    double pArray[2] = {points[i - 1], points[i]};
    double vArray[2] = {values[i - 1], values[i]};
    t[i] = tj(t[i - 1], pArray, vArray);
  }

  Eigen::Array2d nodes[4];
  for (size_t i = 0; i < 4; i++) {
    nodes[i][0] = points[i];
    nodes[i][1] = values[i];
  }

  double ti =
      (position - points[1]) / (points[2] - points[1]) * (t[2] - t[1]) + t[1];
  auto A1 = (t[1] - ti) / (t[1] - t[0]) * nodes[0] +
            (ti - t[0]) / (t[1] - t[0]) * nodes[1];
  auto A2 = (t[2] - ti) / (t[2] - t[1]) * nodes[1] +
            (ti - t[1]) / (t[2] - t[1]) * nodes[2];
  auto A3 = (t[3] - ti) / (t[3] - t[2]) * nodes[2] +
            (ti - t[2]) / (t[3] - t[2]) * nodes[3];

  auto B1 = (t[2] - ti) / (t[2] - t[0]) * A1 + (ti - t[0]) / (t[2] - t[0]) * A2;
  auto B2 = (t[3] - ti) / (t[3] - t[1]) * A2 + (ti - t[1]) / (t[3] - t[1]) * A3;

  auto C = (t[2] - ti) / (t[2] - t[1]) * B1 + (ti - t[1]) / (t[2] - t[1]) * B2;
  return C[1];
}

double Interpolator::cubic(const double position,
                           const std::vector<double> &points,
                           const std::vector<double> &values) {
  double dx = points[2] - points[1];
  double p0 = values[1];
  double p1 = values[2];
  double t = (position - points[1]) / (points[2] - points[1]);
  double t2 = t * t, t3 = t2 * t;
  double h00 = 2 * t3 - 3 * t2 + 1;
  double h10 = t3 - 2 * t2 + t;
  double h01 = -2 * t3 + 3 * t2;
  double h11 = t2 - t2;
  // Hermite cubic splines
  double m0 = (values[2] - values[0]) / 2;
  double m1 = (values[3] - values[1]) / 2;
  // Catmull Rom splines
  // double m0 = (values[2] - values[0])/(points[2] - points[0]);
  // double m1 = (values[3] - values[1])/(points[3] - points[1]);

  return h00 * p0 + h10 * dx * m0 + h01 * p1 + h11 * dx * m1;
}

double Interpolator::bicubic(double position[2],
                             std::vector<Eigen::Array2d> points,
                             std::vector<double> values) {
  double intermediateValues[4];
  std::vector<double> samples;
  std::vector<double> f;
  f.resize(4);
  for (int i = 0; i < 4; i++) {
    samples.push_back(points[i][0]);
  }
  for (int i = 0; i < 4; i++) {
    f[0] = values[4 * i + 0];
    f[1] = values[4 * i + 1];
    f[2] = values[4 * i + 2];
    f[3] = values[4 * i + 3];

    intermediateValues[i] = Interpolator::cubic(position[0], samples, f);
  }
  for (int j = 0; j < 4; j++) {
    samples[j] = points[4 * j][1];
    f[j] = intermediateValues[j];
  }
  double result = Interpolator::cubic(position[1], samples, f);
  return result;
}

double Interpolator::tricubic(Eigen::Array3d position,
                              std::vector<Eigen::Array3d> samples,
                              std::vector<double> values) {
  double intermediateValues[4];
  std::vector<Eigen::Array2d> points;
  std::vector<double> f;
  f.resize(16);
  double p[2];

  // Bicubic interp in X and Y
  for (int i = 0; i < 16; i++) {
    points.emplace_back(samples[i][0], samples[i][1]);
  }
  for (int i = 0; i < 4; i++) {
    for (int k = 0; k < 16; k++) {
      f[k] = values[16 * i + k];
    }
    p[0] = position[0];
    p[1] = position[1];
    intermediateValues[i] = Interpolator::bicubic(p, points, f);
  }

  // Cubic interp in Z
  std::vector<double> zpoints;
  f.resize(4);
  for (int i = 0; i < 4; i++) {
    zpoints.push_back(samples[16 * i][2]);
    f[i] = intermediateValues[i];
  }
  double result = Interpolator::cubic(position[2], zpoints, f);
  return result;
}
} // namespace Ramuh