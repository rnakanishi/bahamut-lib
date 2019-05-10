#include <structures/triangle_mesh.h>
#include <glm/glm.hpp>
#include <iostream>

namespace Ramuh {
TriangleMesh::TriangleMesh() {}
TriangleMesh::~TriangleMesh() {}

uint TriangleMesh::addVertex(glm::vec3 vertex) {
  // if (_vMap.find(vertex) == _vMap.end()) {
  _vertices.emplace_back(vertex);
  return _vertices.size() - 1;
  // _vMap[vertex] = _vertices.size() - 1;
  // }
  // return _vMap[vertex];
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

} // namespace Ramuh