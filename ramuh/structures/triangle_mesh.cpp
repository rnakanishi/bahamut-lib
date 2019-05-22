#include <glm/glm.hpp>
#include <iostream>
#include <structures/triangle_mesh.h>

namespace Ramuh {
TriangleMesh::TriangleMesh() {}
TriangleMesh::~TriangleMesh() {}

uint TriangleMesh::addVertex(glm::vec3 vertex) {
  _vertices.emplace_back(vertex);
  return _vertices.size() - 1;
}

uint TriangleMesh::addFace(glm::ivec3 face) {
  _faces.emplace_back(face);
  return _faces.size() - 1;
}

glm::vec3 TriangleMesh::getVertex(uint index) { return _vertices[index]; }

glm::ivec3 TriangleMesh::getFace(uint index) { return _faces[index]; }

uint TriangleMesh::getVerticesSize() { return _vertices.size(); }

uint TriangleMesh::getFacesSize() { return _faces.size(); }

void remesh();

glm::vec3 TriangleMesh::getCentroid() { return _centroid; }

glm::vec3 TriangleMesh::getBBoxSize() { return _bboxMax - _bboxMin; }

glm::vec3 TriangleMesh::getBboxCenter() { return (_bboxMax + _bboxMin) / 2.0f; }

void TriangleMesh::_computeCentroid() {
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

} // namespace Ramuh