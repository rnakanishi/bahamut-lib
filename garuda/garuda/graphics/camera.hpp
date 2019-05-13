#ifndef __GARUDA_CAMERA_HPP
#define __GARUDA_CAMERA_HPP

#include <glm/vec3.hpp>

namespace Garuda {
class Camera {
public:
  Camera();

  glm::vec3 getPosition();

protected:
  glm::vec3 position;
  glm::vec3 lookAt, lookUp;
  float far, near;
};
} // namespace Garuda

#endif // !__GARUDA_CAMERA_HPP