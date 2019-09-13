#ifndef __LEVIATHAN_PARTICLE_LEVELSET_H__
#define __LEVIATHAN_PARTICLE_LEVELSET_H__

#include <structures/particle_system2.h>
#include <fluids/levelset_fluid2.h>

namespace Leviathan {

class ParticleLevelSet2 : public Ramuh::ParticleSystem2, LevelSetFluid2 {
public:
  ParticleLevelSet2();

protected:
  int _radiusDataId;
};

} // namespace Leviathan

#endif