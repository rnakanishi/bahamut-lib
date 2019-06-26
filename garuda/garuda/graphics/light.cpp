#include <graphics/light.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <sstream>
#include <GLFW/glfw3.h>

namespace Garuda {
Light::Light() {
  _position = glm::vec3(2.f, 0.f, 0.f);
  _color = glm::vec4(1.f, 1.f, 1.f, 1.f);
  _transform = glm::mat4(1.f);
  // _transform = glm::translate(_transform, glm::vec3(-3.f));
}

Light::Light(glm::vec3 position, glm::vec4 color) {
  _position = position;
  _color = color;
  _transform = glm::mat4(1.f);
}

void Light::initialize() {
  unsigned int vbo[2];

  // Create a point to represent the light
  glGenVertexArrays(1, &_vao);
  glBindVertexArray(_vao);
  glGenBuffers(2, vbo);

  glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3), &_position, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

glm::vec3 Light::getPosition() { return _position; }

glm::vec4 Light::getColor() { return _color; }

void Light::setPosition(glm::vec3 position) { _position = position; }

void Light::setPosition(glm::vec4 color) { _color = color; }

glm::mat4 &Light::transformMatrix() { return _transform; }

void Light::draw(Shader shader) {
  shader.useShader();
  unsigned int modelUniform = glGetUniformLocation(shader.getId(), "model");
  glUniformMatrix4fv(modelUniform, 1, GL_FALSE, glm::value_ptr(_transform));

  unsigned int lightUniform = glGetUniformLocation(shader.getId(), "vColor");
  glUniform4f(lightUniform, _color[0], _color[1], _color[2], _color[3]);
  glBindVertexArray(_vao);
  glDrawArrays(GL_POINTS, 0, 1);
}

void Light::illuminate(Shader shader) { illuminate(shader, 0); }

void Light::illuminate(Shader shader, int lightId) {
  std::stringstream str_lightPosition;
  std::stringstream str_lightColor;
  std::stringstream str_transform;

  str_lightPosition << "lights[" << lightId << "].position";
  str_lightColor << "lights[" << lightId << "].color";
  str_transform << "lightTransform[" << lightId << "]";

  int lightUniform =
      glGetUniformLocation(shader.getId(), str_lightPosition.str().c_str());
  int modelUniform =
      glGetUniformLocation(shader.getId(), str_transform.str().c_str());
  glUniform3f(lightUniform, _position[0], _position[1], _position[2]);
  glUniformMatrix4fv(modelUniform, 1, GL_FALSE, glm::value_ptr(_transform));
  lightUniform =
      glGetUniformLocation(shader.getId(), str_lightColor.str().c_str());
  glUniform4f(lightUniform, _color[0], _color[1], _color[2], _color[3]);
}

void Light::move() {
  glm::mat4 transform = glm::mat4(1.f);
  // transform = glm::translate(-_position);
  transform =
      glm::rotate(transform, (float)glfwGetTime(), glm::vec3(0.f, 1.f, 0.f));
  transform = glm::translate(transform, _position);
  _transform = transform;
}
} // namespace Garuda