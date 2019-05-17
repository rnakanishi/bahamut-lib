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

  /***************************************
   * @brief Move camera along its own axis
   *
   * @param direction Direction where camera should move
   *
   * This function moves the camera along the its own axis. The user is unaware
   * about camera positioning in world space, meaning that moving camera
   * forward, the viewing window will move forward.
   ***************************************/
  void moveCamera(glm::vec3 direction);

  /**
   * @brief Rotate camera by pitch and yaw angles
   *
   * @param pitch Rotation angle around x axis
   * @param yaw Rotation angle around y axis
   */
  void rotateCamera(double pitch, double yaw);

  void setLookUp(glm::vec3 lookUp);

  void setLookAt(glm::vec3 lookAt);

  void orthogonalProjection();

  void perspectiveProjection();

  glm::vec3 cameraDirection();

  glm::mat4 viewMatrix();

  glm::mat4 projectionMatrix();

protected:
  glm::vec3 _position;
  glm::vec3 _front, _lookUp, _right, _lookAt;
  glm::mat4 _projection;
  float _far, _near;
  float _pitch, _yaw;
};
} // namespace Garuda

#endif // !__GARUDA_CAMERA_HPP