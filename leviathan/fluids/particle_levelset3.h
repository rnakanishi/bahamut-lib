#ifndef __LEVIATHAN_PARTICLE_LEVELSET_H__
#define __LEVIATHAN_PARTICLE_LEVELSET_H__

#include <structures/particle_system3.h>
#include <fluids/levelset_fluid3.h>
#include <set>

namespace Leviathan {

class ParticleLevelSet3 : public Ramuh::ParticleSystem3, public LevelSetFluid3 {
public:
  ParticleLevelSet3();

  ParticleLevelSet3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain);

  /**
   * @brief Compute Euler advection scheme to the particles. This method does
   * not set velocities values to the particles, so this step should be called
   * beforehand by interpolateVelocityToParticles() method.
   *
   */
  void advectParticles();

  void advectEuler();

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

  void reseedParticles();

  int findCellIdByCoordinate(Eigen::Array3d position);

  void adjustParticleRadius();

  bool correctLevelSetWithParticles();

protected:
  bool _hasEscaped(int pid);

  std::vector<int> _findSurfaceCells(int distacenToSurface);
  std::vector<int> _findSurfaceCells();

  void _seedCells(std::vector<int> &toSeed, std::vector<int> &n);
  void _seedCells(std::vector<int> &toSeed);
  void _seedCells(std::set<int> &toSeed);

  int _maxParticles;

  // Arrays ids
  int _particleRadiusId, _particleVelocityId;
  int _particleSignalId;
  int _particleLevelSetId;

  std::vector<int> _escapedParticles;
};

} // namespace Leviathan

#endif