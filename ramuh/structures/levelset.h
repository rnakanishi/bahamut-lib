#ifndef __RAMUH_LEVELSET_H__
#define __RAMUH_LEVELSET_H__

#include <structures/grid.h>
#include <utils/material.h>

namespace Ramuh {

class LevelSet : public RegularGrid {
public:
  LevelSet();

  ///
  /// Add an implicit region for the fluid
  void addImplicitFunction();

  ///
  /// As long as velocity field is defined over cell faces (MAC grid),
  /// interpolate thos values to cell vertices
  void interpolateVelocitiesToVertices();

  ///
  /// Advect level set according to \f$\phi_t + u\cdot\nabla\phi = 0\f$. This
  /// method assumes that a velocity field is already defined over cell corners
  void integrateLevelSet();

  ///
  /// Define a value for each vertex of the grid correspnoding to the isocontour
  /// of a sphere given its center and radius
  /// \param center center of the sphere
  /// \param radius radius of the sphere
  void addSphereSurface(Vector3d center, double radius);

  void checkCellMaterial();

  void printVertexVelocity();

  void setResolution(Vector3i newResolution) override;

  std::vector<std::vector<double>> &operator[](const int i);

  Vector3d operator()(int i, int j, int k);

protected:
  std::vector<std::vector<std::vector<Vector3d>>> _gradPhi,
      _velocity; // level set gradient and velocity on the corners
  std::vector<std::vector<std::vector<double>>>
      _phi; // Level set stored on the grid corners and its gradient values
};
} // namespace Ramuh

#endif