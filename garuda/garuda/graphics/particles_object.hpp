#ifndef __GARUDA_PARTICLES_OBJECT_HPP__
#define __GARUDA_PARTICLES_OBJECT_HPP__
#include <memory>
#include <graphics/scene_object.hpp>
#include <structures/particle_system3.h>

namespace Garuda {

class ParticlesObject3 : public SceneObject {
public:
  ParticlesObject3(std::shared_ptr<Ramuh::ParticleSystem3> particles);

  void draw(Shader shader) override;

protected:
  std::shared_ptr<Ramuh::ParticleSystem3> _particles;
};

} // namespace Garuda

#endif
