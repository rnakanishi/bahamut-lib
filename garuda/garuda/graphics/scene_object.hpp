#ifndef __GARUDA_SCENE_OBJECT_HPP__
#define __GARUDA_SCENE_OBJECT_HPP__
#include <glm/mat4x4.hpp>
#include <vector>
#include <shader/shader.hpp>

namespace Garuda {

class SceneObject {
public:
  SceneObject();

  /**
   * @brief Creates vbo, vao and ebo buffers.
   *
   **/
  virtual void initialize();

  /**
   * @brief Bind corresponding buffers for rendering
   *
   **/
  virtual void draw(Shader shader) = 0;

protected:
  unsigned int _vbo[3], _vao, _ebo, _tex;
  std::vector<glm::mat4> _modelMatrix, _normalMatrix;
  int _instanceCount;
};

} // namespace Garuda

#endif
