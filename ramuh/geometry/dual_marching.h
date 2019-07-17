#ifndef __RAMUH_DUAL_MARCHING_H__
#define __RAMUH_DUAL_MARCHING_H__
#include <map>
#include <utility>
#include <Eigen/Dense>
#include <geometry/bounding_box.h>

namespace Ramuh {
class DualMarching2 {
public:
  DualMarching2();

  DualMarching2(int resolution);

  /**
   * @brief evaluate the 2d square containing levelset Hermite data. Points
   *should be sorted starting from bottom left and going in a counter clockwise
   *direction. Same rule is applied to edge normals.
   *
   * @param pointIndices vertices indices to be evaluated. This value is used to
   *store the points in the class instance, so the surface can be reconstructed
   *later using reconstruct() function.
   * @param normalLocation normal position for intersecting edges
   * @param normals themselves
   * @return Eigen::Array2d least square position of the new point
   **/
  Eigen::Array2d evaluateSquare(std::pair<int, int> pointIndices,
                                std::vector<Eigen::Array2d> normalLocation,
                                std::vector<Eigen::Vector2d> normals);

  void reconstruct();

private:
  std::map<std::pair<int, int>, Eigen::Array2d> _points;
  int _resolution;
};

class DualMarching3 {
public:
  DualMarching3();

  DualMarching3(Eigen::Array3i resolution);

  Eigen::Array3d evaluateCube(std::tuple<int, int, int> pointIndices,
                              std::vector<Eigen::Array3d> normalLocation,
                              std::vector<Eigen::Vector3d> normals);


  Eigen::Array3d evaluateCube(std::tuple<int, int, int> pointIndices,
                              std::vector<Eigen::Array3d> normalLocation,
                              std::vector<Eigen::Vector3d> normals, 
                              BoundingBox3 cubeLimits);

  void reconstruct();

private:
  std::map<std::tuple<int, int, int>, int> _idMap;
  std::map<int, Eigen::Array3d> _points;
  std::map<int, Eigen::Vector3d> _normals;
  Eigen::Array3i _resolution;
};

} // namespace Ramuh

#endif