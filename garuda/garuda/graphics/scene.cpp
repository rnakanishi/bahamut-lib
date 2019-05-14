#include <graphics/scene.hpp>

namespace Garuda {
Scene::Scene() {
  _cameras.emplace_back(Camera());
  _objects.emplace_back(MeshObject());
}

void Scene::load() {
  _shader.loadVertexShader("./assets/shaders/camera.vert");
  _shader.loadFragmentShader("./assets/shaders/texture.frag");

  for (auto &object : _objects) {
    object.initialize();
    object.loadObjMesh("./assets/3d_models/newdog.obj");
    object.loadTexture();
    object.centerizeObject();
  }
}

void Scene::draw() {
  _shader.useShader();
  for (auto &object : _objects)
    object.draw(_shader);
}

} // namespace Garuda