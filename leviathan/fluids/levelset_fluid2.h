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
   **/
  void advect();

  void advectWeno();

  /**
   * @brief Perform Cubic Interpolated Propagation (CIP) advection for the
   * levelset. The predicted (t+1) gradient is computed using shifted gradients
   * method, and then phi(t+1) is computed. The gradient is corrected using the
   * updated phi value for the first 10 iterations and then after every 50 time
   * steps.
   *
   * Reference paper: "Cubic interpolated pseudo particle method for solving
   * hyperbolic type equations", H. Takawaki, A. Nishiguchi, T. Yabe. JCP 1985
   *
   */
  void advectCip();

  void computeCellsGradient();

  void redistance();

  bool advanceTime();

  void applyCfl();

  virtual void print();

  void trackSurface();

  double interpolateCellScalarData(int dataId, Eigen::Array2d position);

  Eigen::Array2d interpolateCellArrayData(int dataId, Eigen::Array2d position);

protected:
  double __interpolateVelocityU(Eigen::Array2d position);
  double __interpolateVelocityV(Eigen::Array2d position);
  double __interpolateVelocityU(Eigen::Array2d position, double &min,
                                double &max);
  double __interpolateVelocityV(Eigen::Array2d position, double &min,
                                double &max);

  double _dt, _ellapsedDt, _originalDt;
  size_t _velocityId, _phiId, _gradientId;
  bool _isPressure2nd;
  double _tolerance;
  std::vector<int> _surfaceCells;
};

} // namespace Leviathan

#endif