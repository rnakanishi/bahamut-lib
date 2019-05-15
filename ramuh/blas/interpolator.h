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

  static double bilinear();
  static double trilinear(Eigen::Array3d position,
                          std::vector<Eigen::Array3d> points,
                          std::vector<double> values);
};
} // namespace Ramuh

#endif