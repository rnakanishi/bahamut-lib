#ifndef __LEVIATHAN_LEVELSET_FLUID_H__
#define __LEVIATHAN_LEVELSET_FLUID_H__

#include <structures/mac_grid2.h>

namespace Leviathan {

class LevelSetFluid2 : Ramuh::MacGrid2 {
public:
  LevelSetFluid2();

  /*
   * @brief Advect the values int the levelset using the grid velocity. An order
   *for the advection can be given. If no order is set, then semi Lagrangean
   *advection is used (first order accuracy)
   *
   * @param [opt] order of advection accuracy. Default: 1
   **/
  void advect(int order);
  void advect();

protected:
  double _dt;
  bool _isPressure2nd;
  double _maxVelocity[2];
};

} // namespace Leviathan

#endif