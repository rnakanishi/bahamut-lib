#include <structures/line_mesh.hpp>

namespace Ramuh {
LineMesh::LineMesh() {}

int LineMesh::addVertex(Eigen::Array2d position) {
  _verticesPosition.emplace_back(position);
  return _verticesPosition.size() - 1;
}

std::vector<int> LineMesh::addVertices(std::vector<Eigen::Array2d> positions) {
  std::vector<int> indices;
  for (auto position : positions) {
    indices.emplace_back(_verticesPosition.size());
    _verticesPosition.emplace_back(position);
  }
  return indices;
}

int LineMesh::connectVertices(Eigen::Array2i segment) {
  _segments.emplace_back(segment);
  return _segments.size() - 1;
}

std::vector<int>
LineMesh::connectVertices(std::vector<Eigen::Array2i> segments) {
  std::vector<int> indices;
  for (auto segment : segments) {
    indices.emplace_back(_segments.size());
    _segments.emplace_back(segment);
  }
  return indices;
}

int LineMesh::addSegment(Eigen::Array2d origin, Eigen::Array2d ending) {
  int originId = addVertex(origin);
  int endingId = addVertex(ending);
  return connectVertices(Eigen::Array2i(originId, endingId));
}

Eigen::Array2d LineMesh::getVertexPosition(int vertexId) {
  return _verticesPosition[vertexId];
}

Eigen::Array2i LineMesh::getSegmentVertices(int segmentId) {
  return _segments[segmentId];
}

int LineMesh::getVerticesCount() { return _verticesPosition.size(); }

int LineMesh::getSegmentsCount() { return _segments.size(); }

std::vector<Eigen::Array2i> &LineMesh::getSegmentsList() { return _segments; }

std::vector<Eigen::Array2d> &LineMesh::getVerticesList() {
  return _verticesPosition;
}

} // namespace Ramuh
