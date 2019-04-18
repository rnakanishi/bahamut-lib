#include <graphics/triangle_mesh.hpp>

namespace Garuda {
TriangleMesh::TriangleMesh() {}

void TriangleMesh::initialize() {
  glGenVertexArrays(1, &_vao);
  glGenBuffers(1, &_vbo);
}

void TriangleMesh::loadObjMesh() {
  glm::ivec3 face;
  face[0] = addVertex(glm::vec3(-0.5, -0.286, 0.0));
  face[1] = addVertex(glm::vec3(0.5, 0 - 0.286, 0.0));
  face[2] = addVertex(glm::vec3(0.0, 0.574, 0.0));
  addFace(face);

  glBindVertexArray(_vao);
  glBindBuffer(GL_ARRAY_BUFFER, _vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * getVerticesSize(),
               _vertices.data(), GL_STATIC_DRAW);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void TriangleMesh::draw() {
  glBindVertexArray(_vao);
  glDrawArrays(GL_TRIANGLES, 0, 3);
}
} // namespace Garuda