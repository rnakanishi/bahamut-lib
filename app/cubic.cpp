#include <blas/interpolator.h>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <iostream>
#include <iomanip>

double _f(Eigen::Array3d p) {
  return (std::cos(p[0]) + std::sin(p[1]) - std::cos(p[2]));
}
double _f(Eigen::Array2d p) { return (std::cos(p[0]) + std::sin(p[1])); }
double _f(double p) { return std::cos(5 * p); }

void test3D() {
  std::vector<Eigen::Array3d> points;
  std::vector<double> values;
  double h = 1.0 / 40;

  Eigen::Array3d mostBottom(0.1 - h, 0.1 - h, 0.1 - h);
  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        points.emplace_back(mostBottom[0] + i * h, mostBottom[1] + j * h,
                            mostBottom[2] + k * h);
      }
    }
  }

  for (auto p : points)
    values.emplace_back(_f(p));

  Eigen::Array3d target;
  for (size_t i = 0; i < 10; i++) {
    target[0] = (double)(rand() % 100) / 100 * h + 0.1;
    target[1] = (double)(rand() % 100) / 100 * h + 0.1;
    target[2] = (double)(rand() % 100) / 100 * h + 0.1;
    double value = Ramuh::Interpolator::tricubic(target, points, values);

    std::cout << std::setw(8) << target[0] << ", ";
    std::cout << std::setw(8) << target[1] << ", ";
    std::cout << std::setw(8) << target[2] << ":\t\t";
    std::cout << std::setw(8) << value << std::setw(3) << " " << _f(target);
    std::cout << "\t Error: ";
    std::cout << std::fabs(value - _f(target)) << std::endl;
  }
}

void test2D() {
  std::vector<Eigen::Array2d> points;
  std::vector<double> values;
  double h = 1.0 / 40;

  Eigen::Array2d mostBottom(0.1 - h, 0.1 - h);
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
      points.emplace_back(mostBottom[0] + i * h, mostBottom[1] + j * h);
    }
  }

  for (auto p : points)
    values.emplace_back(_f(p));

  double target[2];
  for (size_t i = 0; i < 10; i++) {
    target[0] = (double)(rand() % 1000) / 1000.0 * h + 0.1;
    target[1] = (double)(rand() % 1000) / 1000.0 * h + 0.1;
    // target[0] = (double)(rand() % 100) / 100.0 * 2 * h + 0.1 - h / 2.;
    // target[1] = (double)(rand() % 100) / 100.0 * 2 * h + 0.1 - h / 2.;
    double value = Ramuh::Interpolator::bicubic(target, points, values);

    std::cout << std::setw(8) << target[0] << ", ";
    std::cout << std::setw(8) << target[1] << ":\t\t";
    std::cout << std::setw(8) << std::setprecision(5) << value << std::setw(8)
              << "\t" << _f(Eigen::Array2d(target[0], target[1]));
    std::cout << "\t Error: ";
    std::cout << std::fabs(value - _f(Eigen::Array2d(target[0], target[1])))
              << std::endl;
  }
}

void test1D() {
  std::vector<double> points;
  std::vector<double> values;
  double h = 1.0 / 40;

  points.emplace_back(0.1 - h);
  points.emplace_back(0.1);
  points.emplace_back(0.1 + h);
  points.emplace_back(0.1 + 2 * h);

  for (auto p : points)
    values.emplace_back(_f(p));

  double target;
  for (size_t i = 0; i < 10; i++) {
    target = (double)(rand() % 1000) / 1000.0 * h + 0.1;
    double value = Ramuh::Interpolator::cubic(target, points, values);

    std::cout << std::setw(8) << target << ":\t\t";
    std::cout << std::setw(8) << std::setprecision(5) << value << std::setw(8)
              << ' ' << _f(target);
    std::cout << "\t Error: ";
    std::cout << std::fabs(value - _f(target)) << std::endl;
  }
}

int main(int argc, char const *argv[]) {
  test3D();
  return 0;
}
