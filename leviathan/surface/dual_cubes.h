#ifndef __LEVIATHAN_DUAL_CUBES_H__
#define __LEVIATHAN_DUAL_CUBES_H__

#include <fluids/levelset_fluid3.h>
#include <geometry/dual_marching.h>
#include <string>

namespace Leviathan {

class DualCubes : public LevelSetFluid3 {
public:
  DualCubes(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain);

  DualCubes(LevelSetFluid3 levelset);

  void extractSurface();

  void computeIntersection();

  void computeIntersectionAndNormals();
  Eigen::Array3d computeIntersection(int cell1, int cell2);

  void setFolder(std::string folder);

  void resetFileCounter();
  void resetFileCounter(int value);

  void swapLevelSet(LevelSetFluid3 levelset);

protected:
  bool hasSignChange(double valueA, double valueB);

  int _faceSurfaceNormalId, _faceSurfacePositionId;
  std::string _baseFolder;
  int _count;
  Ramuh::DualMarching3 surface;
};

} // namespace Leviathan

#endif