#ifndef __CARBUNCLE_DUAL_CUBES_H__
#define __CARBUNCLE_DUAL_CUBES_H__
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <geometry/bounding_box.h>
#include <structures/mac_grid3.h>
#include <structures/mac_grid2.h>
#include <fluids/levelset_fluid3.h>
#include <fstream>
#include <sstream>

namespace Carbuncle {

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

class DualCubes3 : public Leviathan::LevelSetFluid3 {
public:
  enum class ParametricSurface : int { SPHERE, CUBE, TORUS, DOUBLE_TORUS };

  /**
   * @brief Construct a new Dual Cubes 3 object. This object has size and
   *resolution predefined as [0,1] and 32^3 respectively.
   *
   **/
  DualCubes3();

  DualCubes3(Eigen::Array3i resolution);

  DualCubes3(Eigen::Array3i resolution, Ramuh::BoundingBox3 domain);

  void initialize(Eigen::Array3d center, double radius);

  void initialize(Eigen::Array3d center, double radius,
                  ParametricSurface surface);

  void analyticNormals(Eigen::Array3d center, double radius,
                       ParametricSurface surface);

  void computeNormals();

  void defineVelocity();

  void computeIntersection();

  void extractSurface();

  void print() {
    auto &phi = getScalarData(_phiId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/weno/3d/" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);

    for (size_t i = 0; i < cellCount(); i++) {
      file << phi[i] << "\n";
    }
    file.close();
  }

private:
  bool signChange(double valueA, double valueB);

protected:
  std::vector<std::vector<std::vector<double>>> _vertices;
  std::vector<std::vector<std::vector<Eigen::Vector3d>>> _ufaceNormals,
      _vfaceNormals, _wfaceNormals;
  std::vector<std::vector<std::vector<Eigen::Array3d>>> _ufaceLocation,
      _vfaceLocation, _wfaceLocation;
  Eigen::Array3d _h;
  Eigen::Array3i _resolution;
};
} // namespace Carbuncle

#endif
