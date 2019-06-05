#ifndef __GARUDA_LIGHT_HPP
#define __GARUDA_LIGHT_HPP

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

namespace Garuda {
class Light {
public:
  Light();

  glm::vec3 getPosition();

  glm::vec4 getColor();

  void setPosition();

  void setPosition();

protected:
  glm::vec3 _position;
  glm::vec4 _color;
};
} // namespace Garuda

#endif // !__GARUDA_LIGHT_HPP