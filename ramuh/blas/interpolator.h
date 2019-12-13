#ifndef __RAMUH_BLAS_INTERPOLATOR_H
#define __RAMUH_BLAS_INTERPOLATOR_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <geometry/vector2.h>

namespace Ramuh {
class Interpolator {
public:
  Interpolator();

  /**
   * @brief Performs a linear interpolation of the values along the axis.
   *
   * @param position target for the interpolation
   * @param points point coordinates surrounding the value. This verification is
   * NOT performed, so it is assumed user gives correct values
   * @param values values of the above given coordinates
   * @return double interpolated value
   */
  static double linear(double &position, std::vector<double> &points,
                       std::vector<double> &values);

  /**
   * @brief Performs a bilinear interpolation of the values. For each
   * coordinate, a separate interpolation is done. First for all pair of values
   * along x-coordinate. The these values are used to perform interpolation on
   * y-axis.
   *
   * Data should be given in this order: X and then Y
   *
   * @param position target for the interpolation
   * @param points point coordinates of the square surrounding the value. This
   * verification is NOT performed, so it is assumed user gives correct values
   * @param values values of the above given coordinates
   * @return double interpolated value
   */
  static double bilinear(double position[2],
                         std::vector<Eigen::Array2d> &points,
                         std::vector<double> &values);

  /**
   * @brief Performs a trilinear interpolation of the values. For each
   * coordinate, a separate interpolation is done. First for all pair of values
   * along x-coordinate. The these values are used to perform interpolation on
   * y-axis. And finally the resulting values are used to interpolate along
   * z-axis.
   *
   * Data should be given in this order: X, Y and then Z
   *
   * @param position target for the interpolation
   * @param points point coordinates of the cube surrounding the value. This
   * verification is NOT performed, so it is assumed user gives correct values
   * @param values values of the above given coordinates
   * @return double interpolated value
   */
  static double trilinear(Eigen::Array3d &position,
                          std::vector<Eigen::Array3d> &points,
                          std::vector<double> &values);

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

  static double bicubic(double position[2], std::vector<Eigen::Array2d> &points,
                        std::vector<double> values);

  static double tricubic(Eigen::Array3d position,
                         std::vector<Eigen::Array3d> &points,
                         std::vector<double> values);

  static double shepard(double position, std::vector<double> samples,
                        std::vector<double> values);

  /**
   * @brief Computes shepard interpolation (aka inverse distance interpolation)
   * for points in 3D space. This method assumes that points and their values
   * are passed in the same order.
   * If a sample point is too close from the query point, the that value is
   * returned without further computation
   *
   * @param position query point for the interpolation
   * @param samples all points that will be used to interpolate
   * @param values the respective values to be interpolated
   * @return double final interpolated value at queried position
   */
  static double shepard(Eigen::Array3d position,
                        std::vector<Eigen::Array3d> &samples,
                        std::vector<double> values);

  static Eigen::Vector3d shepard(Eigen::Array3d position,
                                 std::vector<Eigen::Array3d> &samples,
                                 std::vector<Eigen::Vector3d> vectors);

  static double shepard(Eigen::Array2d position,
                        std::vector<Eigen::Array2d> &samples,
                        std::vector<double> values);

  static Eigen::Vector2d shepard(Eigen::Array2d position,
                                 std::vector<Eigen::Array2d> &samples,
                                 std::vector<Eigen::Vector2d> vectors);

  static double rbf(Eigen::Array3d position,
                    std::vector<Eigen::Array3d> &samples,
                    std::vector<double> values);

  static double catmullRom(double position, std::vector<double> points,
                           std::vector<double> values);

  static double camullRom2(Eigen::Array2d position,
                           std::vector<Eigen::Array2d> &points,
                           std::vector<double> values);

  static double camullRom2(Eigen::Array3d position,
                           std::vector<Eigen::Array3d> &points,
                           std::vector<double> values);

  static double closestPoint(Eigen::Array2d position,
                             std::vector<Eigen::Array2d> &points,
                             std::vector<double> values);
  static Eigen::Vector2d closestPoint(Eigen::Array2d position,
                                      std::vector<Eigen::Array2d> &points,
                                      std::vector<Eigen::Vector2d> values);
};
} // namespace Ramuh

#endif