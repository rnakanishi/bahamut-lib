#include <graphics/camera.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/vec4.hpp>

#include <iostream>
namespace Garuda {

Camera::Camera() {
  _position = glm::vec3(0.0f, 0.0f, 3.0f);
  _lookUp = glm::vec3(0.0f, 1.0f, 0.0f);
  _lookAt = glm::vec3(0.f);
  _front = glm::normalize(_lookAt - _position);

  _projection = glm::mat4();
  _near = 0.1f;
  _far = 100.0f;

  _yaw = -90.f;
  _pitch = 0.f;

  _right = glm::normalize(glm::cross(_lookUp, _position));
  setLookUp(glm::cross(_position, _right));
}

glm::vec3 Camera::getPosition() { return _position; }

void Camera::setPosition(glm::vec3 position) { _position = position; }

void Camera::moveCamera(glm::vec3 direction) {
  // Transform back to original position
  // TODO: Align axis to origin to perform general movementation
  float speed = 0.05f;
  glm::mat4 translate = glm::translate(glm::mat4(1.f), -_position);
  glm::mat4 rotate = viewMatrix();
  glm::mat4 alignment = glm::transpose(rotate * translate);

  glm::vec4 position = glm::vec4(_position[0], _position[1], _position[2], 1.f);
  position = glm::inverse(alignment) * position;
  position += glm::vec4(direction, 1.f) * speed;
  position = alignment * position;
  _position = glm::vec3(position);
}

void Camera::setLookAt(glm::vec3 lookAt) {
  _lookAt = lookAt;
  _front = glm::normalize(lookAt - _position);
  _right = glm::normalize(glm::cross(_front, _lookUp));
  setLookUp(glm::cross(_right, _front));
}

void Camera::rotateCamera(double pitch, double yaw) {
  glm::mat4 translate = glm::translate(glm::mat4(1.f), -_position);
  glm::mat4 rotate = viewMatrix();
  glm::mat4 alignment = glm::transpose(rotate * translate);

  _pitch += pitch;
  _yaw += yaw;
  _pitch = std::min(std::max(_pitch, -89.f), 89.f);

  glm::vec4 front = glm::vec4(_front, 1.f);
  front = glm::inverse(alignment) * front;
  rotate = glm::mat4();
  rotate = glm::rotate(rotate, (float)pitch, glm::vec3(1.0f, 0.f, 0.f));
  rotate = glm::rotate(rotate, (float)yaw, glm::vec3(0.0f, 1.f, 0.f));
  front = glm::transpose(rotate) * front;
  front = alignment * front;
  _front = glm::normalize(glm::vec3(front.x, front.y, front.z));
}

void Camera::setLookUp(glm::vec3 lookUp) { _lookUp = glm::normalize(lookUp); }

void Camera::orthogonalProjection() {
  _projection = glm::ortho(0.0f, 600.0f, 0.0f, 600.0f, _near, _far);
}

void Camera::perspectiveProjection() {
  _projection = glm::perspective(glm::radians(45.f), 4.f / 4.f, _near, _far);
}

glm::vec3 Camera::cameraDirection() { return glm::normalize(_front); }

glm::mat4 Camera::viewMatrix() {
  return glm::lookAt(_position, _position + _front, _lookUp);
}

glm::mat4 Camera::projectionMatrix() { return _projection; }

} // namespace Garuda