#ifndef __LEVIATHAN_LEVELSET_FLUID_H__
#define __LEVIATHAN_LEVELSET_FLUID_H__

#include <structures/mac_grid2.h>

namespace Leviathan {

class LevelSetFluid2 : public Ramuh::MacGrid2 {
public:
  LevelSetFluid2();

  LevelSetFluid2(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain);

  /*
   * @brief Advect the values int the levelset using the grid velocity. An order
   *for the advection can be given. If no order is set, then semi Lagrangean
   *advection is used (first order accuracy)
   *
   * @param [opt] order of advection accuracy. Default: 1
   **/
  void advect();

  void advectWeno();

  void computeCellsGradient();

  void redistance();

  virtual void print();

protected:
  double __interpolateVelocityU(Eigen::Array3d position);
  double __interpolateVelocityV(Eigen::Array3d position);
  double __interpolateVelocityU(Eigen::Array3d position, double &min,
                                double &max);
  double __interpolateVelocityV(Eigen::Array3d position, double &min,
                                double &max);

  double _dt;
  size_t _velocityId, _phiId;
  bool _isPressure2nd;
  double _tolerance;
};

} // namespace Leviathan

#endif