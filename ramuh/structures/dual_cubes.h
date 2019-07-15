#ifndef __RAMUH_DUAL_CUBES_H__
#define __RAMUH_DUAL_CUBES_H__
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <geometry/bounding_box.h>

namespace Ramuh {

class DualCubes2 {
public:
  DualCubes2();

  DualCubes2(int size);

  void initialize(Eigen::Array2d center, double radius);

  void extractSurface();

  void printCells();

private:
  bool signChange(double valueA, double valueB);

protected:
  std::vector<std::vector<double>> _cells, _vertices;
  std::vector<std::vector<Eigen::Vector2d>> _ufaceNormals, _vfaceNormals;
  std::vector<std::vector<Eigen::Array2d>> _ufaceLocation, _vfaceLocation;
  double _h;
  double _resolution;
};

class DualCubes3 {
public:
  /**
   * @brief Construct a new Dual Cubes 3 object. This object has size and
   *resolution predefined as [0,1] and 32^3 respectively.
   *
   **/
  DualCubes3();

  DualCubes3(Eigen::Array3i resolution);

  DualCubes3(Eigen::Array3i resolution, BoundingBox3 domain);

  void initialize(Eigen::Array3d center, double radius);

  void extractSurface();

  void printCells();

private:
  bool signChange(double valueA, double valueB);

protected:
  std::vector<std::vector<std::vector<double>>> _cells, _vertices;
  std::vector<std::vector<std::vector<Eigen::Vector3d>>> _ufaceNormals,
      _vfaceNormals, _wfaceNormals;
  std::vector<std::vector<std::vector<Eigen::Array3d>>> _ufaceLocation,
      _vfaceLocation, _wfaceLocation;
  Eigen::Array3d _h;
  BoundingBox3 _domain;
  Eigen::Array3i _resolution;
};
} // namespace Ramuh

#endif