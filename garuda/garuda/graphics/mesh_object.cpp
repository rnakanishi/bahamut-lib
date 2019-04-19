#include <graphics/mesh_object.hpp>
#include <utils/mesh_reader.hpp>

namespace Garuda {
MeshObject::MeshObject() {}

void MeshObject::initialize() {
  glGenVertexArrays(1, &_vao);
  glGenBuffers(2, _vbo);
  glGenBuffers(1, &_ebo);
}

void MeshObject::loadObjMesh(const char *objPath) {
  // glm::ivec3 face, color;
  // face[0] = addVertex(glm::vec3(-0.5, -0.286, 0.0));
  // face[1] = addVertex(glm::vec3(0.5, -0.286, 0.0));
  // face[2] = addVertex(glm::vec3(0.0, 0.574, 0.0));
  // addFace(face);
  // face[2] = addVertex(glm::vec3(0.0, -0.860, 0.0));
  // addFace(face);
  _vertexColor.emplace_back(glm::vec3(1.0, 0.0, 0.0));
  _vertexColor.emplace_back(glm::vec3(0.0, 1.0, 0.0));
  _vertexColor.emplace_back(glm::vec3(0.0, 0.0, 1.0));
  _vertexColor.emplace_back(glm::vec3(0.0, 1.0, 1.0));

  MeshReader::readObj(*this, objPath);

  glBindVertexArray(_vao);

  // Assign vertex position buffer
  glBindBuffer(GL_ARRAY_BUFFER, _vbo[0]);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * getVerticesSize(),
               _vertices.data(), GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
  glEnableVertexAttribArray(0);

  // Assign vertices indices to triangles
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(glm::ivec3) * getFacesSize(),
               _faces.data(), GL_STATIC_DRAW);

  // Assign color to each vertex
  glBindBuffer(GL_ARRAY_BUFFER, _vbo[1]);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * getVerticesSize(),
               _vertexColor.data(), GL_STATIC_DRAW);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
  glEnableVertexAttribArray(1);

  // Ensure no buffer is binded
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void MeshObject::draw(Shader shader) {
  shader.useShader();
  glBindVertexArray(_vao);
  // glDrawArrays(GL_TRIANGLES, 0, 3 * getVerticesSize());
  glDrawElements(GL_TRIANGLES, 3 * getFacesSize(), GL_UNSIGNED_INT, 0);
}

} // namespace Garuda