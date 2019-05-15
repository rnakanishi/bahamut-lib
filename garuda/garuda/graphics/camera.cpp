#include <graphics/camera.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <iostream>
namespace Garuda {

Camera::Camera() {
  _position = glm::vec3(0.0f, 1.0f, 3.0f);
  _lookUp = glm::vec3(0.0f, 1.0f, 0.0f);
  _lookAt = glm::vec3(0.0f, 0.0f, 0.0f);

  _projection = glm::mat4();
  _near = 0.1f;
  _far = 100.0f;

  _right = glm::cross(_lookUp, _position);
}

glm::vec3 Camera::getPosition() { return _position; }

void Camera::setPosition(glm::vec3 position) { _position = position; }

void Camera::setLookAt(glm::vec3 lookAt) { _lookAt = lookAt; }

void Camera::setLookUp(glm::vec3 lookUp) { _lookUp = lookUp; }

void Camera::orthogonalProjection() {
  _projection = glm::ortho(0.0f, 600.0f, 0.0f, 600.0f, _near, _far);
}

void Camera::perspectiveProjection() {
  _projection = glm::perspective(glm::radians(45.f), 4.f / 4.f, _near, _far);
}

glm::vec3 Camera::cameraDirection() {
  return glm::normalize(_lookAt - _position);
}

glm::mat4 Camera::viewMatrix() {
  return glm::lookAt(_position, _lookAt, _lookUp);
}

glm::mat4 Camera::projectionMatrix() { return _projection; }

} // namespace Garuda