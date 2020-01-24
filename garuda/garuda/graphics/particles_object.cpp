#include <graphics/particles_object.hpp>

namespace Garuda {
ParticlesObject3::ParticlesObject3(
    std::shared_ptr<Ramuh::ParticleSystem3> particles) {
  _particles.swap(particles);
}

} // namespace Garuda
