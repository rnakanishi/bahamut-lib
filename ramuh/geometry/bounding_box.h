#ifndef __RAMUH_BOUNDING_BOX_H__
#define __RAMUH_BOUNDING_BOX_H__
#include <Eigen/Dense>

namespace Ramuh {

class BoundingBox3 {
public:
  BoundingBox3();

  BoundingBox3(Eigen::Array3d min, Eigen::Array3d max);

  BoundingBox3(double min, double max);

  /**
   * @brief creates a new instance of a bounding box with min coordinates at
   *origin and unit side.
   *
   * @return BoundingBox3
   **/
  static BoundingBox3 unitBox();

  /**
   * @brief Set the lowest corner coordinate for the bounding box
   *
   * @param new min coordinate
   */
  void setMin(Eigen::Array3d min);

  /**
   * @brief Set the highest corner coordinate for the bounding box
   *
   * @param new max coordinate
   */
  void setMax(Eigen::Array3d max);

  /**
   * @brief Returns this bpunding box total size with double precision. This
   *value is computed taking the difference between each corner coordinate.
   *
   * @return Eigen::Array3d Each coordinate of the array corresponds to the
   *respective side's size.
   **/
  Eigen::Array3d getSize();

  /**
   * @brief Get the Min coodinate of the bounding box
   *
   * @return Eigen::Array3d coordinate of the lowest corner
   */
  Eigen::Array3d getMin();

  /**
   * @brief Get the Max coordinate of the bounding box
   *
   * @return Eigen::Array3d  coordinate of the highest corner
   */
  Eigen::Array3d getMax();

  /**
   * @brief Get the Center coordinate of the bounding box. This coordinate is
   * computed by taking the average coordinates of the bounding box corners
   *
   * @return Eigen::Array3d center coordinate of the bounding box
   */
  Eigen::Array3d getCenter();

  /**
   * @brief Clamp the point passed as parameter to the bounding box. If any
   * coordinate goes out the box, it is set to the box limits insteads
   *
   * @param point to be clamped
   * @return Eigen::Array3d clamped point
   */
  Eigen::Array3d clamp(Eigen::Array3d point);

  /**
   * @brief Check either the bounding box contains the point given as parameter.
   * If the point is inside the box, then true is returned
   *
   * @return true if the point is inside the box
   * @return false otherwise
   */
  bool contains(Eigen::Array3d);

private:
  Eigen::Array3d _min, _max;
};

class BoundingBox2 {
public:
  BoundingBox2();

  BoundingBox2(Eigen::Array2d min, Eigen::Array2d max);

  BoundingBox2(double min, double max);

  /**
   * @brief creates a new instance of a bounding box with min coordinates at
   *origin and unit side.
   *
   * @return BoundingBox2
   **/
  static BoundingBox2 unitBox();

  /**
   * @brief Returns this bpunding box total size with double precision
   *
   * @return Eigen::Array2d Each coordinate of the array corresponds to the
   *respective side's size.
   **/
  Eigen::Array2d size();

  Eigen::Array2d min();

  Eigen::Array2d max();

  Eigen::Array2d center();

  void setMin(Eigen::Array2d newMin);

  void setMax(Eigen::Array2d newMax);

  Eigen::Array2d clamp(Eigen::Array2d point);

  bool contains(Eigen::Array2d);

  bool contains(BoundingBox2 box);

private:
  Eigen::Array2d _min, _max;
};

class BoundingBox1 {
public:
  BoundingBox1();

  BoundingBox1(double min, double max);

  /**
   * @brief creates a new instance of a bounding box with min coordinates at
   *origin and unit side.
   *
   * @return BoundingBox2
   **/
  static BoundingBox1 unitBox();

  /**
   * @brief Returns this bpunding box total size with double precision
   *
   * @return double Each coordinate of the array corresponds to the
   *respective side's size.
   **/
  double size();

  double min();

  double max();

  double center();

  double clamp(double point);

  bool contains(double);

private:
  double _min, _max;
};

} // namespace Ramuh

#endif
