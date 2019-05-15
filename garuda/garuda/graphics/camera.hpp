#ifndef __GARUDA_CAMERA_HPP
#define __GARUDA_CAMERA_HPP

#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>

namespace Garuda {
class Camera {
public:
  Camera();

  glm::vec3 getPosition();

  void setPosition(glm::vec3 position);

  void translate(glm::vec3 direction);

  void setLookUp(glm::vec3 lookUp);

  void setLookAt(glm::vec3 lookAt);

  void orthogonalProjection();

  void perspectiveProjection();

  glm::vec3 cameraDirection();

  glm::mat4 viewMatrix();

  glm::mat4 projectionMatrix();

protected:
  glm::vec3 _position;
  glm::vec3 _lookAt, _lookUp, _right;
  glm::mat4 _projection;
  float _far, _near;
};
} // namespace Garuda

#endif // !__GARUDA_CAMERA_HPP