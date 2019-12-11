#ifndef __LEVIATHAN_DUAL_SQUARES_H__
#define __LEVIATHAN_DUAL_SQUARES_H__

#include <fluids/levelset_fluid2.h>
#include <structures/mac_grid2.h>
#include <geometry/dual_marching.h>
#include <string>

namespace Leviathan {
class DualSquares : public LevelSetFluid2 {
public:
  DualSquares(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain);

  DualSquares(LevelSetFluid2 levelset);

  void swapLevelSet(LevelSetFluid2 &levelset);

  void extractSurface();

  void computeIntersectionAndNormals();

  void setFolder(std::string folder);

  void resetFileCounter();

  void resetFileCounter(int value);

  bool hasSignChange(double valueA, double valueB);

protected:
  int _faceSurfaceNormalId, _faceSurfacePositionId;
  std::string _baseFolder;
  int _count;
  Ramuh::DualMarching2 surface;
};

} // namespace Leviathan

#endif
