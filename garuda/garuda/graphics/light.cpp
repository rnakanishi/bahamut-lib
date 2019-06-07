#include <graphics/light.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <GLFW/glfw3.h>

namespace Garuda {
Light::Light() {
  _position = glm::vec3(5.f, 5.f, 5.f);
  _color = glm::vec4(1.f, 1.f, 1.f, 1.f);
  _transform = glm::mat4(1.f);
  // _transform = glm::translate(_transform, glm::vec3(-3.f));
}

void Light::initialize() {
  unsigned int vbo;
  // Create a point to represent the light
  glGenVertexArrays(1, &_vao);
  glGenBuffers(1, &vbo);
  glBindVertexArray(_vao);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
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

  glBindVertexArray(_vao);
  glDrawArrays(GL_POINTS, 0, 1);
}

void Light::illuminate(Shader shader) {
  int lightUniform = glGetUniformLocation(shader.getId(), "light");
  int modelUniform = glGetUniformLocation(shader.getId(), "lightTransform");

  glUniform3f(lightUniform, _position[0], _position[1], _position[2]);
  glUniformMatrix4fv(modelUniform, 1, GL_FALSE, glm::value_ptr(_transform));
}

void Light::move() {
  _transform = glm::rotate(_transform, 0.05f, glm::vec3(-1.f, 1.f, -1.f));
}
} // namespace Garuda