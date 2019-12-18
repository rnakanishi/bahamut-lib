#include <structures/line_mesh.hpp>

namespace Ramuh {
LineMesh::LineMesh() {}

int LineMesh::addVertex(Eigen::Array2d position) {
  int vId = _verticesPosition.size();
  _verticesPosition.emplace_back(position);
  _vertSegments[vId] = std::vector<int>();
  return vId;
}

std::vector<int> LineMesh::addVertices(std::vector<Eigen::Array2d> positions) {
  std::vector<int> indices;
  for (auto position : positions) {
    indices.emplace_back(addVertex(position));
  }
  return indices;
}

int LineMesh::connectVertices(int vertex1, int vertex2) {
  return connectVertices(Eigen::Array2i(vertex1, vertex2));
}

int LineMesh::connectVertices(Eigen::Array2i segment) {
  int segId = hasConnection(segment[0], segment[1]);
  if (segId < 0) {
    segId = _segments.size();
    _segments.emplace_back(segment);
    _vertSegments[segment[0]].emplace_back(segId);
    _vertSegments[segment[1]].emplace_back(segId);
  }
  return segId;
}

std::vector<int>
LineMesh::connectVertices(std::vector<Eigen::Array2i> segments) {
  std::vector<int> indices;
  for (auto segment : segments) {
    indices.emplace_back(connectVertices(segment));
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

std::vector<int> &LineMesh::getVertexSegments(int vertexId) {
  return _vertSegments[vertexId];
}

int LineMesh::hasConnection(int vertex1, int vertex2) {
  auto connections = _vertSegments[vertex1];
  for (auto segmentId : connections) {
    auto segment = _segments[segmentId];
    if (segment[0] == vertex2 || segment[1] == vertex2)
      return segmentId;
  }
  connections = _vertSegments[vertex2];
  for (auto segmentId : connections) {
    auto segment = _segments[segmentId];
    if (segment[0] == vertex1 || segment[1] == vertex1)
      return segmentId;
  }
  return -1;
}

int LineMesh::getNumberOfConnections(int vertexId) {
  return _vertSegments[vertexId].size();
}

std::vector<int> LineMesh::getAdjacentVertices(int vertexId) {
  std::vector<int> neighbors;
  auto segIds = getVertexSegments(vertexId);

  for (auto segId : segIds) {
    auto vertIds = getSegmentVertices(segId);
    if (vertIds[0] != vertexId)
      neighbors.emplace_back(vertIds[0]);
    else
      neighbors.emplace_back(vertIds[1]);
  }
  return neighbors;
}

void LineMesh::disconnectVertices(int vertex1, int vertex2) {
  // Find which segment is between those two vertices
  auto segId = hasConnection(vertex1, vertex2);
}

} // namespace Ramuh
