#ifndef __LEVIATHAN_PARTICLE_LEVELSET_H__
#define __LEVIATHAN_PARTICLE_LEVELSET_H__

#include <structures/particle_system2.h>
#include <fluids/levelset_fluid2.h>

namespace Leviathan {

class ParticleLevelSet2 : public Ramuh::ParticleSystem2, public LevelSetFluid2 {
public:
  ParticleLevelSet2();

  ParticleLevelSet2(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain);

  /**
   * @brief Compute Euler advection scheme to the particles. This method does
   * not set velocities values to the particles, so this step should be called
   * beforehand by interpolateVelocityToParticles() method.
   *
   */
  void advectParticles();

  /**
   * @brief Interpolate grid face velocities to the particles using bilinear
   * interpolation method.
   *
   */
  void interpolateVelocityToParticles();

  /**
   * @brief From the tracked surface cells, compute where particles should be
   * seeded. Compute bfs distance from the interface and cells within 3
   * manhattan distance from surface receive particles.
   * Note that particles receive random position and signal.
   *
   */
  void seedParticlesNearSurface();

  /**
   * @brief Procede with the attraction phase for particles, so their signal
   * should match levelset's signal. A max iteration of 15 is performed and if
   * the particle ends up on the wrong side, it is deleted.
   *
   */
  void attractParticles();

  void adjustParticleRadius();

  void correctLevelSetWithParticles();

protected:
  bool _hasEscaped(int pid);

  int _particleRadiusId, _particleVelocityId;
  int _particleSignalId;
  int _particleLevelSetId;

  std::vector<int> _escapedParticles;
};

} // namespace Leviathan

#endif