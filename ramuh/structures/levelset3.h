#ifndef __RAMUH_LEVELSET3_H__
#define __RAMUH_LEVELSET3_H__

#include <structures/grid3.h>
#include <utils/material.h>

namespace Ramuh {

class LevelSet3 : public RegularGrid3 {
public:
  LevelSet3();

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
  /// TODO: Change to semi lagrangean method
  void integrateLevelSet();

  ///
  /// Define a value for each vertex of the grid correspnoding to the isocontour
  /// of a sphere given its center and radius
  /// \param center center of the sphere
  /// \param radius radius of the sphere
  void addSphereSurface(Vector3d center, double radius);
  void addCubeSurface(Vector3d lower, Vector3d upper);

  void checkCellMaterial();

  void printVertexVelocity();

  void setResolution(Vector3i newResolution) override;

  std::vector<std::vector<double>> &operator[](const int i);

  Vector3d operator()(int i, int j, int k);

protected:
  // TODO: Change to Matrix3 type
  std::vector<std::vector<std::vector<Vector3d>>> _gradPhi,
      _velocity; // level set gradient and velocity on the corners
  std::vector<std::vector<std::vector<double>>>
      _phi; // Level set stored on the grid corners and its gradient values
};
} // namespace Ramuh

#endif