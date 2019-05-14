#ifndef __GARUDA_SCENE_HPP
#define __GARUDA_SCENE_HPP

#include <graphics/mesh_object.hpp>

namespace Garuda {
class Scene {
public:
  Scene();

protected:
  std::vector<MeshObject> _objects;
  std::vector<Camera> _cameras;
};
} // namespace Garuda

#endif