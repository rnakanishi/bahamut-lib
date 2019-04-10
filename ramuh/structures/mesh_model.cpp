#include <structures/mesh_model.h>

namespace Ramuh {
MeshModel3::MeshModel3() {}
MeshModel3::~MeshModel3() {}

uint MeshModel3::addVertex(glm::vec3 vertex) {
  _vertices.emplace_back(vertex);
  return _vertices.size() - 1;
}

uint MeshModel3::addFace(glm::ivec3 face) {
  _faces.emplace_back(face);
  return _faces.size() - 1;
}

glm::vec3 MeshModel3::getVertex(uint index) { return _vertices[index]; }

glm::ivec3 MeshModel3::getFace(uint index) { return _faces[index]; }

uint MeshModel3::getVerticesSize() { return _vertices.size(); }

uint MeshModel3::getFacesSize() { return _faces.size(); }

} // namespace Ramuh