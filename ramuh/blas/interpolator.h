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

  /**
   * @brief Cubic Hermite Splines interpolation for arbitrary interval
   *  Assumes that only four points are passed as parameter. It computes cubic
   *splines using conventional tangents. TODO: Monotone interpolation is yet to
   *be implemented
   *
   * @param position points between second and third sample point
   * @param points Vector with fout points
   * @param values Function values for given points
   * @return double Interpolated value at position
   **/
  static double cubic(const double position, const std::vector<double> &points,
                      const std::vector<double> &values);

  static double bicubic(double position[2], std::vector<Eigen::Array2d> points,
                        std::vector<double> values);

  static double tricubic(Eigen::Array3d position,
                         std::vector<Eigen::Array3d> points,
                         std::vector<double> values);

  static double shepard(double position, std::vector<double> samples,
                        std::vector<double> values);

  static double shepard(Eigen::Array3d position,
                        std::vector<Eigen::Array3d> samples,
                        std::vector<double> values);

static double rbf(Eigen::Array3d position,
                        std::vector<Eigen::Array3d> samples,
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