#ifndef __RAMUH_BLAS_INTERPOLATOR_H
#define __RAMUH_BLAS_INTERPOLATOR_H

#include <vector>
#include <Eigen/Dense>

namespace Ramuh {
class Interpolator {
public:
  Interpolator();

  static double linear(double position, std::vector<double> points,
                       std::vector<double> values);

  static double bilinear(double position[2], std::vector<Eigen::Array2d> points,
                         std::vector<double> values);

  static double trilinear(Eigen::Array3d position,
                          std::vector<Eigen::Array3d> points,
                          std::vector<double> values);

  static double cubic(double position, std::vector<double> points,
                      std::vector<double> values);

  static double catmullRom(double position, std::vector<double> points,
                           std::vector<double> values);

  static double camullRom2(Eigen::Array2d position,
                           std::vector<Eigen::Array2d> points,
                           std::vector<double> values);

  static double camullRom2(Eigen::Array3d position,
                           std::vector<Eigen::Array3d> points,
                           std::vector<double> values);
};
} // namespace Ramuh

#endif