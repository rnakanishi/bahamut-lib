#include <graphics/light.hpp>

namespace Garuda {
Light::Light() {
  _position = glm::vec3(5.f, 5.f, 5.f);
  _color = glm::vec4(1.f, 1.f, 1.f, 1.f);
}

} // namespace Garuda