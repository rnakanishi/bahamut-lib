#ifndef __CARBUNCLE_NORMAL_PARTICLES_HPP__
#define __CARBUNCLE_NORMAL_PARTICLES_HPP__

#include <structures/particle_system3.h>
#include <fluids/levelset_fluid3.h>
#include <Eigen/Dense>
#include <map>

namespace Carbuncle {

class NormalParticles3 : public Ramuh::ParticleSystem3 {
public:
  NormalParticles3();

  int seedParticlesOverSurface(Leviathan::LevelSetFluid3 levelset);

  void print();

protected:
  std::map<int, int> normalPair;
};

} // namespace Carbuncle

#endif
