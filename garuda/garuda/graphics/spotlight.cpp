#include <graphics/light.hpp>

namespace Garuda {
SpotLight::SpotLight() : Light() {
  _angle = 45.f;
  _direction = glm::vec3();
}

void SpotLight::setDirection(glm::vec3 direction) { _direction = direction; }

void SpotLight::setAngle(float angle) { _angle = angle; }

glm::vec3 SpotLight::getDirection() { return _direction; }

float SpotLight::getAngle() { return _angle; }

} // namespace Garuda