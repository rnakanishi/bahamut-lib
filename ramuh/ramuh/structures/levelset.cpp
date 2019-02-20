#include <structures/levelset.h>
#include <utils/macros.h>

namespace Ramuh {

LevelSet::LevelSet() : RegularGrid() {
  _phi.resize((_resolution.x() + 1) * (_resolution.y() + 1) *
              (_resolution.z() + 1));
}

void LevelSet::addImplicitFunction() { NOT_IMPLEMENTED(); }

void LevelSet::interpolateVelocitiesToVertices() { NOT_IMPLEMENTED(); }

void LevelSet::integrateLevelSet() { NOT_IMPLEMENTED(); }

void LevelSet::solvePressure() {

  // Compute velocity divergent over cell center

  // Solve pressure Poisson equation

  // Correct velocity through pressure gradient
}
} // namespace Ramuh