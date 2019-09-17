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
   * @brief Returns this bpunding box total size with double precision
   *
   * @return Eigen::Array3d Each coordinate of the array corresponds to the
   *respective side's size.
   **/
  Eigen::Array3d size();

  Eigen::Array3d min();

  Eigen::Array3d max();

  Eigen::Array3d center();

  Eigen::Array3d clamp(Eigen::Array3d point);

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
