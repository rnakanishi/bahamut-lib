#ifndef __GARUDA_TRIANGLE_MESH_HPP__
#define __GARUDA_TRIANGLE_MESH_HPP__

#include <glad/glad.h>
#include <structures/mesh_model.h>

namespace Garuda {
class TriangleMesh : public Ramuh::MeshModel3 {
public:
  TriangleMesh();

  void initialize();

  void draw();

  void loadObjMesh();

protected:
  unsigned int _vbo, _vao, _ebo;
};

} // namespace Garuda

#endif