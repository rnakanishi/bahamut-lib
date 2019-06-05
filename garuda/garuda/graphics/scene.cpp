#include <graphics/scene.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>

namespace Garuda {
Scene::Scene() {
  _cameras.emplace_back(Camera());
  _objects.emplace_back(MeshObject());
  activeCamera = 0;
  _ambientLight = glm::vec3(0.6, 0.8, 1.0);
}

void Scene::load() {
  _shader.loadVertexShader("./assets/shaders/light.vert");
  _shader.loadFragmentShader("./assets/shaders/light.frag");

  for (auto &object : _objects) {
    object.initialize();
    object.loadObjMesh("./assets/3d_models/cow2.obj");
    object.loadTexture();
    object.centerizeObject();
  }

  _cameras[0].perspectiveProjection();
}

void Scene::draw() {
  _shader.useShader();
  Camera activeCamera = _cameras[0];

  // glClearColor(_ambientLight[0], _ambientLight[1], _ambientLight[2], 1.0f);
  glClearColor(0.2f, 0.2f, 0.2f, 1.f);
  int ambientLigheLocation =
      glGetUniformLocation(_shader.getId(), "ambientLight");
  glUniform3f(ambientLigheLocation, _ambientLight[0], _ambientLight[1],
              _ambientLight[2]);

  int projectionUniform = glGetUniformLocation(_shader.getId(), "projection");
  glUniformMatrix4fv(projectionUniform, 1, GL_FALSE,
                     glm::value_ptr(activeCamera.projectionMatrix()));

  int viewUniform = glGetUniformLocation(_shader.getId(), "view");
  glUniformMatrix4fv(viewUniform, 1, GL_FALSE,
                     glm::value_ptr(activeCamera.viewMatrix()));

  for (auto &object : _objects)
    object.draw(_shader);

  for (auto &camera : _cameras)
    camera.draw();
}

Camera &Scene::getActiveCamera() { return _cameras[activeCamera]; }

MeshObject &Scene::getObject(int index) { return _objects[index]; }

} // namespace Garuda