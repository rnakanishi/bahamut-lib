#ifndef __LEVIATHAN_LEVELSET_FLUID_2_H__
#define __LEVIATHAN_LEVELSET_FLUID_2_H__

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
  void advectUpwind();

  void advectSemiLagrangian();

  void advectRungeKutta3();

  void computeCentralGradient();

  void computeWenoGradient();

  void advectWeno();

  void wenoAdvection();

  /**
   * @brief Perform Cubic Interpolated Propagation (CIP) advection for the
   * levelset. The predicted (t+1) gradient is computed using shifted
   * gradients method, and then phi(t+1) is computed. The gradient is
   * corrected using the updated phi value for the first 10 iterations and
   * then after every 50 time steps.
   *
   * Reference paper: "Cubic interpolated pseudo particle method for solving
   * hyperbolic type equations", H. Takawaki, A. Nishiguchi, T. Yabe. JCP 1985
   *
   */
  void advectCip();

  void computeCellsGradient();

  void computeCellVelocity();

  void redistanceSimple();

  void redistance();

  /**
   * @brief Evaluates if timestep moved forward, accorgindly to cfl condition.
   * If all cfl pieces arec ompleted, then return true.
   *
   * @return true If time step evolved
   * @return false If more substeps are needed
   */
  bool advanceTime();

  void applyCfl();

  virtual void print();

  // TODO: remove this fucntion from particlelevelset2
  virtual int findCellIdByCoordinate(Eigen::Array2d position);

  std::vector<int> findSurfaceCells();
  std::vector<int> findSurfaceCells(double surfaceDistance);

  bool isSurfaceCell(int cellId);

  std::vector<int> findCellNeighbors(int cellId, int distance = 1);

  std::vector<int> trackSurface();

  void merge(LevelSetFluid2 &levelset);

protected:
  double __interpolateVelocityU(Eigen::Array2d position);
  double __interpolateVelocityV(Eigen::Array2d position);
  double __interpolateVelocityU(Eigen::Array2d position, double &min,
                                double &max);
  double __interpolateVelocityV(Eigen::Array2d position, double &min,
                                double &max);

  double _dt, _ellapsedDt, _originalDt;
  size_t _faceVelocityId, _phiId, _gradientId;
  size_t _cellVelocityId;
  bool _isPressure2nd;
  double _tolerance;
  int _surfaceCellCount;
  std::vector<bool> _isSurfaceCell;
  std::vector<int> _surfaceCellIds;
};

} // namespace Leviathan

#endif