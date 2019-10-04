#ifndef __LEVIATHAN_PARTICLE_LEVELSET_H__
#define __LEVIATHAN_PARTICLE_LEVELSET_H__

#include <structures/particle_system2.h>
#include <fluids/levelset_fluid2.h>
#include <set>

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

  int findCellIdByCoordinate(Eigen::Array2d position);

  void adjustParticleRadius();

  bool correctLevelSetWithParticles();

  void sortParticles();

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
  std::map<int, std::vector<int>> _particlesInCell;
};

} // namespace Leviathan

#endif