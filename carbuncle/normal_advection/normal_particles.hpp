#ifndef __CARBUNCLE_NORMAL_PARTICLES_HPP__
#define __CARBUNCLE_NORMAL_PARTICLES_HPP__

#include <structures/particle_system3.h>
#include <fluids/levelset_fluid3.h>
#include <fluids/particle_levelset3.h>
#include <Eigen/Dense>
#include <map>

namespace Carbuncle {

class NormalParticles3 : public Leviathan::ParticleLevelSet3 {
public:
  NormalParticles3();

  int seedParticlesOverSurface(Leviathan::LevelSetFluid3 levelset);

  std::map<int, int> &getPairMap();

  void estimateCellNormals(Leviathan::LevelSetFluid3 &levelset);

  void advectParticles(double time);

protected:
  std::map<int, int> normalPair;
  int _particleVelocityId;
};

} // namespace Carbuncle

#endif
