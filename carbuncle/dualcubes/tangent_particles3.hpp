#ifndef __CARBUNCLE_TANGENT_PARTICLES_3_HPP__
#define __CARBUNCLE_TANGENT_PARTICLES_3_HPP__

#include <fluids/particle_levelset3.h>
#include <structures/triangle_mesh.h>
#include <map>

namespace Carbuncle {
class TangentParticles3 : public Leviathan::ParticleLevelSet3 {
public:
  TangentParticles3();

  void seedParticleOverSurface();

  void assignPair(int particleId, int pairId);

  void removeParticle(int pid);

  int getParticlePairId(int pid);

  void estimateCellNormals(Leviathan::LevelSetFluid3 &levelset);

protected:
  std::map<int, int> tangentPair;
};
} // namespace Carbuncle

#endif