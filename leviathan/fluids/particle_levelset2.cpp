#include <fluids/particle_levelset2.h>

namespace Leviathan {

ParticleLevelSet2::ParticleLevelSet2() {
  _radiusDataId = newScalarLabel("radius");
}

} // namespace Leviathan
