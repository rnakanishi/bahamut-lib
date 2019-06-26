#ifndef __GARUDA_SCENE_HPP__
#define __GARUDA_SCENE_HPP__
////
#include <glm/glm.hpp>
#include <graphics/camera.hpp>
#include <graphics/mesh_object.hpp>
#include <graphics/light.hpp>

namespace Garuda {
class Scene {
public:
  Scene();

  void load();

  void draw();

  Camera &getActiveCamera();

  /**
   * @brief Get the MeshObject object
   *
   * @param index
   * @return MeshObject& referred by the index value
   */
  MeshObject &getObject(int index);

protected:
  std::vector<MeshObject> _objects;
  std::vector<Camera> _cameras;
  std::vector<Light> _lights;
  glm::vec3 _ambientLight;
  int _activeCameraId;
  Shader _shader;
  Shader _lampShader;
};
} // namespace Garuda

#endif