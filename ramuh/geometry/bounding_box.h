#ifndef __RAMUH_BOUNDING_BOX_H__
#define __RAMUH_BOUNDING_BOX_H__
#include <Eigen/Dense>

namespace Ramuh {

class BoundingBox2 {
public:
private:
  Eigen::Array2d _min, _max;
};

class BoundingBox3 {
public:
  BoundingBox3();

  BoundingBox3(Eigen::Array3d min, Eigen::Array3d max);

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

private:
  Eigen::Array3d _min, _max;
};

} // namespace Ramuh

#endif
