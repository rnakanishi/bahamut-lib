#ifndef __GARUDA_SCENE_HPP
#define __GARUDA_SCENE_HPP
////
#include <glm/glm.hpp>
#include <graphics/camera.hpp>
#include <graphics/mesh_object.hpp>

namespace Garuda {
class Scene {
public:
  Scene();

  void load();

  void draw();

  Camera &getActiveCamera();

protected:
  std::vector<MeshObject> _objects;
  std::vector<Camera> _cameras;
  int activeCamera;
  Shader _shader;
};
} // namespace Garuda

#endif