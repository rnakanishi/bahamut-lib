#ifndef __RAMUH_LEVELSET_H__
#define __RAMUH_LEVELSET_H__

#include <structures/grid.h>

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
  /// Computes velocities divergent in cell center and then solves pressure
  /// Poisson equation. After, updates the velocity values.
  void solvePressure();

  void setResolution(Vector3i newResolution) override;

  std::vector<std::vector<double>> &operator[](const int i);

  Vector3d operator()(int i, int j, int k);

protected:
  double dt; // Time step
  std::vector<std::vector<std::vector<Vector3d>>>
      _gradPhi; // Level set stored on the grid corners and its gradient values
  std::vector<std::vector<std::vector<double>>>
      _phi; // Level set stored on the grid corners and its gradient values
};
} // namespace Ramuh

#endif