#include <graphics/mesh_object.hpp>
#include <utils/mesh_reader.hpp>
#include <glm/vec2.hpp>
#include <iostream>

namespace Garuda {
MeshObject::MeshObject() {}

void MeshObject::initialize() {
  glGenVertexArrays(1, &_vao);
  glGenBuffers(2, _vbo);
  glGenBuffers(1, &_ebo);
}

void MeshObject::loadObjMesh(const char *objPath) {
  MeshReader::readObj(*this, objPath);
  _computeCentroid();

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

void MeshObject::addTextureCoordinate(glm::vec2 texCoord) {
  _vertexTexture.emplace_back(glm::vec3(texCoord, 0.0));
}

glm::vec3 MeshObject::getCentroid() { return _centroid; }

glm::vec3 MeshObject::getBBoxSize() { return _bboxMax - _bboxMin; }

glm::vec3 MeshObject::getBboxCenter() { return (_bboxMax + _bboxMin) / 2.0f; }

void MeshObject::_computeCentroid() {
  _centroid = glm::vec3(0.0);
  _bboxMax = glm::vec3(-1e8);
  _bboxMin = glm::vec3(1e8);
  for (auto vertex : _vertices) {
    _centroid += vertex;

    if (_bboxMax[0] < vertex[0])
      _bboxMax[0] = vertex[0];
    if (_bboxMax[1] < vertex[1])
      _bboxMax[1] = vertex[1];
    if (_bboxMax[2] < vertex[2])
      _bboxMax[2] = vertex[2];

    if (_bboxMin[0] > vertex[0])
      _bboxMin[0] = vertex[0];
    if (_bboxMin[1] > vertex[1])
      _bboxMin[1] = vertex[1];
    if (_bboxMin[2] > vertex[2])
      _bboxMin[2] = vertex[2];
  }
  _centroid /= (float)_vertices.size();
}

} // namespace Garuda