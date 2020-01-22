#ifndef __CARBUNCLE_TANGENT_PARTICLES2_HPP__
#define __CARBUNCLE_TANGENT_PARTICLES2_HPP__

#include <fluids/particle_levelset2.h>
#include <structures/line_mesh.hpp>
#include <map>

namespace Carbuncle {
class TangentParticles2 : public Leviathan::ParticleLevelSet2 {
public:
  TangentParticles2();

  TangentParticles2(Leviathan::LevelSetFluid2 levelset);

  void defineRotationCenter(Eigen::Array2d center);

  int seedParticlesOverSurface(Leviathan::LevelSetFluid2 levelset);
  int seedParticlesOverSurface(Leviathan::LevelSetFluid2 levelset,
                               Ramuh::LineMesh mesh);

  int seedParticleOverSegment(Eigen::Array2d origin, Eigen::Array2d ending,
                              double cellH, int n = 20);

  void removeParticle(int pid) override;

  std::map<int, int> &getPairMap();

  int getParticlePairId(int pid);

  void assignPair(int particleId, int tangentId);

  void estimateCellNormals(Leviathan::LevelSetFluid2 &levelset);

  void defineParticlesVelocity();

  void fixLevelsetGradients(Leviathan::LevelSetFluid2 &levelset);

  void advectParticles() override;

  void advectParticles(std::function<void(void)> defineParticlesVelocity);

  Ramuh::LineMesh extractSurface(Leviathan::LevelSetFluid2 &levelset,
                                 bool onlyConnections = false);

  static TangentParticles2 mergeParticles(TangentParticles2 p1,
                                          TangentParticles2 p2);

protected:
  std::map<int, int> tangentPair;
  Eigen::Array2d rotationCenter;
};

} // namespace Carbuncle

#endif