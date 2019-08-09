#ifndef __RAMUH_DUAL_MARCHING_H__
#define __RAMUH_DUAL_MARCHING_H__
#include <map>
#include <vector>
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

  /**
   * @brief This method evaluates a 3D cube containing levelset Hermite data.
   *Points indices are used as a key identifier so every vertex is unique.
   *Intersection points of the cube and its normals should be given as a vector
   *with same index for each data. After all cubes of the domain are evalueated,
   *reconstruct() method should be called to generate proper surface polygons.
   *
   * @param pointIndices unique index for each surface point
   * @param normalLocation 3D position of the intersection points
   * @param normals normalLocation correspondent normal vector
   * @return Eigen::Array3d surface point that minimizes QEF function
   **/
  Eigen::Array3d evaluateCube(std::tuple<int, int, int> pointIndices,
                              std::vector<Eigen::Array3d> normalLocation,
                              std::vector<Eigen::Vector3d> normals);

  /**
   * @brief This method evaluates a 3D cube containing levelset Hermite data.
   *Points indices are used as a key identifier so every vertex is unique.
   *Intersection points of the cube and its normals should be given as a vector
   *with same index for each data. A limiting cube could be given so if the new
   *point falls off this bounds, then coordinates are clamped accordingly. After
   *all cubes of the domain are evalueated, reconstruct() method should be
   *called to generate proper surface polygons.
   *
   * @param pointIndices unique index for each surface point
   * @param normalLocation 3D position of the intersection points
   * @param normals normalLocation correspondent normal vector
   * @param cubeLimits boundary limiting cube
   * @return Eigen::Array3d surface point that minimizes QEF function
   **/
  Eigen::Array3d evaluateCube(std::tuple<int, int, int> pointIndices,
                              std::vector<Eigen::Array3d> normalLocation,
                              std::vector<Eigen::Vector3d> normals,
                              BoundingBox3 cubeLimits);

  /**
   * @brief After all cubes of the domain are processed, then this method is
   *called to attach all neighbor cubes together. A different file is written
   *for each call containing the surface vertices position and its normals, as
   *well as the faces for the surface.
   *
   **/
  void reconstruct();

  /**
   * @brief merge the contents of a previously instantiated cube into this
   * instance. This method should be used with caution.
   *
   * @param cube the data to be copied
   */
  void merge(DualMarching3 cube);

  std::vector<Eigen::Array3d> &getPoints();

  std::vector<Eigen::Vector3d> &getNormals();

  std::map<std::tuple<int, int, int>, int> &getIdMap();

private:
  std::map<std::tuple<int, int, int>, int> _idMap;
  std::vector<Eigen::Array3d> _points;
  std::vector<Eigen::Vector3d> _normals;
  Eigen::Array3i _resolution;
};

} // namespace Ramuh

#endif