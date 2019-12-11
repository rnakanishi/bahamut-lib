#ifndef __CARBUNCLE_NORMAL_PARTICLES2_HPP__
#define __CARBUNCLE_NORMAL_PARTICLES2_HPP__

#include <structures/particle_system2.h>
#include <fluids/levelset_fluid2.h>
#include <fluids/particle_levelset2.h>
#include <Eigen/Dense>
#include <map>

namespace Carbuncle {

class NormalParticles2 : public Leviathan::ParticleLevelSet2 {
public:
  NormalParticles2();

  int seedParticlesOverSurface(Leviathan::LevelSetFluid2 levelset);

  std::map<int, int> &getPairMap();

  void estimateCellNormals(Leviathan::LevelSetFluid2 &levelset);

  void defineParticlesVelocity();

  void fixLevelsetGradients(Leviathan::LevelSetFluid2 &levelset);

  void advectParticles() override;

protected:
  std::map<int, int> normalPair;
  int _particleVelocityId;
};

} // namespace Carbuncle

#endif
