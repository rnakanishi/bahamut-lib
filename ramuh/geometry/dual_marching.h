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

  DualMarching2(Eigen::Array2i resolution);

  int convertKey(Eigen::Array2i index);
  int convertKey(int i, int j);
  Eigen::Array2i convertKey(int index);

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
  Eigen::Array2d evaluateSquare(Eigen::Array2i pointIndices,
                                std::vector<Eigen::Array2d> normalLocation,
                                std::vector<Eigen::Vector2d> normals);
  Eigen::Array2d evaluateSquare(Eigen::Array2i pointIndices,
                                std::vector<Eigen::Array2d> normalLocation,
                                std::vector<Eigen::Vector2d> normals,
                                BoundingBox2 squareLimits);

  void reconstruct();
  void reconstruct(std::vector<Eigen::Array2i> connections);

  void merge(DualMarching2 square);

  std::vector<Eigen::Array2d> &getPoints();

  std::vector<Eigen::Vector2d> &getNormals();

  std::map<int, int> &getIdMap();

protected:
  bool _consistentNormals(std::vector<int> ids);

  bool _hasConnection(Eigen::Array2i tuple1, Eigen::Array2i tuple2);

  void _buildConnectionMap(std::vector<Eigen::Array2i> connections);

private:
  std::map<int, int> _idMap;
  std::vector<Eigen::Array2d> _points;
  std::vector<Eigen::Vector2d> _normals;
  std::map<int, std::vector<int>> _connections;
  Eigen::Array2i _resolution;
};

/*************************************************************************
 *************************************************************************
 * @brief DualMarching3 class
 *
 */
class DualMarching3 {
public:
  DualMarching3();

  DualMarching3(Eigen::Array3i resolution);

  int convertKey(Eigen::Array3i index);
  int convertKey(int, int, int);
  Eigen::Array3i convertKey(int index);

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
  Eigen::Array3d evaluateCube(Eigen::Array3i pointIndices,
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
   * @param cellIndex unique index for each surface point
   * @param normalLocation 3D position of the intersection points
   * @param normals normalLocation correspondent normal vector
   * @param cubeLimits boundary limiting cube
   * @return Eigen::Array3d surface point that minimizes QEF function
   **/
  Eigen::Array3d evaluateCube(Eigen::Array3i cellIndex,
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
  void reconstruct(std::vector<std::pair<int, int>> connections);

  /**
   * @brief merge the contents of a previously instantiated cube into this
   * instance. This method should be used with caution.
   *
   * @param cube the data to be copied
   */
  void merge(DualMarching3 cube);

  std::vector<Eigen::Array3d> &getPoints();

  std::vector<Eigen::Vector3d> &getNormals();

  std::map<int, int> &getIdMap();

protected:
  bool _consistentNormals(std::vector<int> ids);

  bool _hasConnection(Eigen::Array3i tuple1, Eigen::Array3i tuple2);

  void _buildConnectionMap(std::vector<std::pair<int, int>> connections);

private:
  std::map<int, int> _idMap;
  std::vector<Eigen::Array3d> _points;
  std::vector<Eigen::Vector3d> _normals;
  Eigen::Array3i _resolution;
  std::map<int, std::vector<int>> _connections;
};

} // namespace Ramuh

#endif