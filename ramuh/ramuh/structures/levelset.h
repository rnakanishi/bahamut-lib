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

protected:
  double dt;                  // Time step
  std::vector<Vector3d> _phi; // Level set stored on the grid corners
};
} // namespace Ramuh