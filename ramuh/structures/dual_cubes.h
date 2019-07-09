#ifndef __RAMUH_DUAL_CUBES_H__
#define __RAMUH_DUAL_CUBES_H__
#include <Eigen/Dense>
#include <vector>

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
  DualCubes3();

  DualCubes3(int size);

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
  double _h;
  double _resolution;
};
} // namespace Ramuh

#endif