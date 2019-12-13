#ifndef __CARBUNCLE_TANGENT_PARTICLES2_HPP__
#define __CARBUNCLE_TANGENT_PARTICLES2_HPP__

#include <fluids/particle_levelset2.h>
#include <map>

namespace Carbuncle {
class TangentParticles2 : public Leviathan::ParticleLevelSet2 {
public:
  TangentParticles2();

  TangentParticles2(Leviathan::LevelSetFluid2 levelset);

  int seedParticlesOverSurface(Leviathan::LevelSetFluid2 levelset);

  std::map<int, int> &getPairMap();

  void estimateCellNormals(Leviathan::LevelSetFluid2 &levelset);

  void defineParticlesVelocity();

  void fixLevelsetGradients(Leviathan::LevelSetFluid2 &levelset);

  void advectParticles() override;

  void extractSurface(Leviathan::LevelSetFluid2 &levelset);

protected:
  std::map<int, int> tangentPair;
};

} // namespace Carbuncle

#endif