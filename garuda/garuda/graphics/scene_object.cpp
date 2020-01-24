#include <graphics/scene_object.hpp>

namespace Garuda {

void SceneObject::initialize() {
  glGenVertexArrays(1, &_vao);
  glGenBuffers(3, _vbo);
  glGenBuffers(1, &_ebo);

  _instanceCount = 1;
  _modelMatrix.emplace_back(glm::mat4(1.0));
  _normalMatrix.emplace_back(glm::mat4(1.0));
}

} // namespace Garuda