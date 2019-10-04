#ifndef __LEVIATHAN_LEVELSET_FLUID_H__
#define __LEVIATHAN_LEVELSET_FLUID_H__

#include <structures/mac_grid3.h>
#include <fstream>
#include <sstream>

namespace Leviathan {

class LevelSetFluid3 : public Ramuh::MacGrid3 {
public:
  LevelSetFluid3();
  LevelSetFluid3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain);

  /**
   * @brief compute semi Lagrangean advection for levelset values: compute
   * inverse time integration and then interpolate value where the cell position
   * would be in previous time.
   *
   * This method assumes cell velocity and cell gradient are previously
   * computed.If not, please call computeCellVelocity() method and
   * computeCellsGradient() method. This method may not perform as expected if
   * those values are not updated.
   */
  void advectSemiLagrangean();

  void advectEuler();

  void advectRungeKutta3();

  void advectMacCormack();

  /**
   * @brief Perform time integration using TVD Runge Kutta algorithm for time
   * discretization and WENO scheme for spatial discretization.
   * Velocities and gradients will be updated in between RK3 substeps.
   *
   * For initial velocity and gradient:
   * Cell velocity and cell gradient are assumed to be computed a priori. If
   * not, please call computeCellVelocity() method and computeCellsGradient()
   * method. This method may not perform as expected if those values are not
   * updated.
   */
  void advectWeno();

  void advectCip();

  /**
   * @brief Performs a time integration using RK1 (Euler) as temporal integrator
   * and upwind schemes for spatial discretization.
   *
   * Cell velocity and cell gradient are assumed to be computed a priori. If
   * not, please call computeCellVelocity() method and computeCellsGradient()
   * method. This method may not perform as expected if those values are not
   * updated.
   */
  void advectUpwind();

  /**
   * @brief Perform upwind/downwind finite differences schemes based on velocity
   * direction. This method modifies values for internal gradient vector data
   *
   * Cell's velocity should be updated using computeCellVelocity() method before
   * this method is called. If not, computed gradient may not be accurate.
   *
   */
  void computeCellsGradient();

  /**
   * @brief Compute WENO scheme for discretizing spatial derivative of the
   * levelset
   *
   */
  void computeWenoGradient();

  /**
   * @brief Using the velocities defined over cell faces, take the average for
   * each coordinate and assign to cell center position
   *
   */
  void computeCellVelocity();

  /**
   * @brief Reinitializes levelset values, so the entire domain is nearer a
   * signed distance field. The algorithm used is based on solving a pde \phi_t
   * + w \cdot \nabla\phi = S(\phi^0),  where w is the gradient vector computed
   * using finite differences scheme.
   *
   * References:
   *   - "A level set method for computing solutions to incompressible two phase
   * flow", M. Sussman, P. Smereka, S. Osher, JCP 1994
   *   - "A remark on computing distance functions", G. Russo and P.
   * Smereka, JCP 2000
   *
   */
  void redistance();

  /**
   * @brief Should be used together with applyCfl() method. After splitting, or
   * not, the timestep into smaller sub steps, verify if all substeps have been
   * computed. If all of them are complete, reset dt to its original state and
   * return true.
   *
   * @return true if all substeps completed.
   * @return false otherwisen
   */
  bool advanceTime();

  /**
   * @brief Computes CFL condition number and, if necessary, splits the dt into
   * accoding pieces to fit CFL number
   *TODO: implement a variation of this method to take cfl number as parameter
   */
  void applyCfl();

  /**
   * @brief Find all cells that are close to the interface. This is made
   * computing the product between two cells. If the resulting value is
   * negative, than the cell is possibly an interface cell. All the cells are
   * internally marked and also their ids are returned in a std::vector
   * container
   *
   * @return std::vector<int> contains all the cell ids that are part of the
   * interface
   */
  std::vector<int> trackSurface();

  std::vector<int> findSurfaceCells(int distacenToSurface);
  std::vector<int> findSurfaceCells();

protected:
  size_t _cellVelocityId, _phiId, _cellGradientId;
  size_t _faceVelocityId;
  bool _isPressure2nd;
  double _dt, _originalDt, _ellapsedDt;
  double _tolerance;

  /// All cells that are part of the interface are marked as true
  std::vector<bool> _isSurfaceCell;
  std::vector<int> _surfaceCellIds;
};

} // namespace Leviathan

#endif