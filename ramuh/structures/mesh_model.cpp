#include <mesh_model.h>

namespace Ramuh {
MeshModel3::MeshModel3() {}
MeshModel3::~MeshModel3() {}

uint MeshModel3::addVertex(Eigen::Vector3d vertex) {
  _vertices.emplace_back(vertex);
  return _vertices.size() - 1;
}

uint MeshModel3::addFace(Eigen::Vector3i face) {
  _faces.emplace_back(face);
  return _faces.size() - 1;
}

} // namespace Ramuh