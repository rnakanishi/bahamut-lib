#ifndef __GARUDA_TRIANGLE_MESH_HPP__
#define __GARUDA_TRIANGLE_MESH_HPP__

#include <glad/glad.h>
#include <structures/triangle_mesh.h>

namespace Garuda {
class TriangleMesh : public Ramuh::TriangleMesh {
public:
  TriangleMesh();

  /**
   * @brief Creates vbo, vao and ebo buffers.
   *
   **/
  void initialize();

  /**
   * @brief Bind corresponding buffers for rendering
   *
   **/
  void draw();

  /**
   * @brief Load an obj mesh and assemble buffers created before.
   *
   **/
  void loadObjMesh(const char *objPath);

protected:
  // TODO: change vbo to accept any size
  unsigned int _vbo[2], _vao, _ebo;
  std::vector<glm::vec3> _vertexColor;
  std::vector<glm::vec3> _vertexTexture;
};

} // namespace Garuda

#endif