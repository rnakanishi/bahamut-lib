#include <graphics/scene.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <gui/event_handler.hpp>

namespace Garuda {
Scene::Scene() {
  _cameras.emplace_back(Camera());
  _objects.emplace_back(MeshObject());
  _lights.emplace_back(Light());
  // _lights.emplace_back(
  // Light(glm::vec3(-1.f, .8f, 1.75f), glm::vec4(0.4f, 0.95f, 0.5f, 1.f)));
  // _lights.emplace_back(
  // Light(glm::vec3(-.2f, -.7f, -0.4f), glm::vec4(.7f, 0.2f, 0.25f, 1.f)));
  _activeCameraId = 0;
  _ambientLight = glm::vec3(0.75, 0.75, 0.75);
}

void Scene::load() {
  _shader.loadVertexShader("./assets/shaders/phong.vert");
  _shader.loadFragmentShader("./assets/shaders/phong.frag");
  _lampShader.loadVertexShader("./assets/shaders/lamp.vert");
  _lampShader.loadFragmentShader("./assets/shaders/lamp.frag");

  for (auto &object : _objects) {
    object.initialize();
    object.loadObjMesh("./assets/3d_models/torus.obj");
    // object.loadObjMesh("./assets/3d_models/bridge/untitled.obj");
    object.loadTexture();
    // object.centerizeObject();
  }

  for (auto light : _lights) {
    light.initialize();
  }

  _cameras[0].perspectiveProjection();
}

void Scene::draw() {
  Camera activeCamera = _cameras[_activeCameraId];

  _shader.useShader();
  // Setting ambient light
  glClearColor(_ambientLight[0], _ambientLight[1], _ambientLight[2], 1.0f);
  int ambientLightLocation =
      glGetUniformLocation(_shader.getId(), "ambientLight");
  glUniform3f(ambientLightLocation, _ambientLight[0], _ambientLight[1],
              _ambientLight[2]);

  // Setting camera properties
  activeCamera.draw(_shader);
  int projectionUniform = glGetUniformLocation(_shader.getId(), "projection");
  glUniformMatrix4fv(projectionUniform, 1, GL_FALSE,
                     glm::value_ptr(activeCamera.projectionMatrix()));
  int viewUniform = glGetUniformLocation(_shader.getId(), "view");
  glUniformMatrix4fv(viewUniform, 1, GL_FALSE,
                     glm::value_ptr(activeCamera.viewMatrix()));
  int camPositionUniform = glGetUniformLocation(_shader.getId(), "camPosition");
  auto camPosition = activeCamera.getPosition();
  glUniform3f(camPositionUniform, camPosition[0], camPosition[1],
              camPosition[2]);

  for (auto &object : _objects)
    object.draw(_shader);

  for (int i = 0; i < _lights.size(); i++) {
    _lights[i].move();
    _lights[i].illuminate(_shader, i);
  }
  unsigned int nLightUniform =
      glGetUniformLocation(_shader.getId(), "n_activeLights");
  glUniform1i(nLightUniform, _lights.size());

  for (auto &camera : _cameras)
    camera.draw(_shader);

  _lampShader.useShader();
  projectionUniform = glGetUniformLocation(_lampShader.getId(), "projection");
  glUniformMatrix4fv(projectionUniform, 1, GL_FALSE,
                     glm::value_ptr(activeCamera.projectionMatrix()));
  viewUniform = glGetUniformLocation(_lampShader.getId(), "view");
  glUniformMatrix4fv(viewUniform, 1, GL_FALSE,
                     glm::value_ptr(activeCamera.viewMatrix()));
  _lampShader.useShader();
  for (auto &light : _lights)
    light.draw(_lampShader);
}

Camera &Scene::getActiveCamera() { return _cameras[_activeCameraId]; }

MeshObject &Scene::getObject(int index) { return _objects[index]; }

} // namespace Garuda