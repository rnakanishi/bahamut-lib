#ifndef __RAMUH_LEVELSET2_H__
#define __RAMUH_LEVELSET2_H__

#include <structures/grid2.h>
#include <utils/material.h>

namespace Ramuh {

class LevelSet2 : public RegularGrid2 {
public:
  LevelSet2();

  ///
  /// Add an implicit region for the fluid
  void addImplicitFunction();

  ///
  /// Advect level set according to \f$\phi_t + u\cdot\nabla\phi = 0\f$. This
  /// method assumes that a velocity field is already defined over cell corners
  /// TODO: Change to semi lagrangean method
  void integrateLevelSet();

  ///
  /// Level set often gets unusual beahvior when advected. To avoid unphysical
  /// behavior, a redistancing function is applied so the level set keep its
  /// property as signed distance
  void redistance();

  ///
  /// Define a value for each vertex of the grid correspnoding to the isocontour
  /// of a sphere given its center and radius
  /// \param center center of the sphere
  /// \param radius radius of the sphere
  void addSphereSurface(Vector2d center, double radius);
  void addCubeSurface(Vector2d lower, Vector2d upper);

  ///
  /// Walks by all level set values (distance field) and check its signal. If
  /// negative, the cell is marked as a fluid cell. This procedure improves
  /// computation of the pressure Poisson system.
  void checkCellMaterial();

  void printVertexVelocity();

  void setResolution(Vector2i newResolution) override;

  std::vector<double> &operator[](const int i);

  Vector2d operator()(int i, int j);

protected:
  // TODO: Change to Matrix3 type
  std::vector<std::vector<Vector2d>> _gradPhi,
      _velocity; // level set gradient and velocity on the corners
  std::vector<std::vector<double>>
      _phi; // Level set stored on the grid corners and its gradient values
};
} // namespace Ramuh

#endif