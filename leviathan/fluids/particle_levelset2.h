#ifndef __LEVIATHAN_PARTICLE_LEVELSET_H__
#define __LEVIATHAN_PARTICLE_LEVELSET_H__

#include <structures/particle_system2.h>
#include <fluids/levelset_fluid2.h>

namespace Leviathan {

class ParticleLevelSet2 : public Ramuh::ParticleSystem2, public LevelSetFluid2 {
public:
  ParticleLevelSet2();

  ParticleLevelSet2(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain);

  void advectParticles();

  void interpolateVelocityToParticles();

protected:
  int _particleRadiusId, _particleVelocityId;
};

} // namespace Leviathan

#endif