#include <structures/line_mesh.hpp>
#include <utils/log.hpp>
#include <string>

namespace Ramuh {
LineMesh::LineMesh() {}

int LineMesh::addVertex(Eigen::Array2d position) {
  return addVertex(position, Eigen::Vector2d(0.0, 0.0));
}

int LineMesh::addVertex(Eigen::Array2d position, Eigen::Vector2d normal) {
  int vId;
  if (_vertexIdQueue.empty()) {
    vId = _verticesPosition.size();
    _verticesPosition.emplace_back(0.);
    _verticesNormals.emplace_back(0., 0.);
    _activeVertices.emplace_back(true);
    _isCorner.emplace_back(false);
  } else {
    vId = _vertexIdQueue.front();
    _vertexIdQueue.pop();
  }
  _activeVertices[vId] = true;
  _verticesPosition[vId] = position;
  _verticesNormals[vId] = normal;
  _vertSegments[vId] = std::vector<int>();
  return vId;
}

std::vector<int> LineMesh::addVertices(std::vector<Eigen::Array2d> positions) {
  std::vector<Eigen::Vector2d> normals(positions.size(),
                                       Eigen::Vector2d(0., 0.));
  return addVertices(positions, normals);
}

std::vector<int> LineMesh::addVertices(std::vector<Eigen::Array2d> positions,
                                       std::vector<Eigen::Vector2d> normals) {
  if (positions.size() != normals.size()) {
    Log::raiseError(
        "LineMesh::AddVertices: Vertices count doesnt match normal count");
    return std::vector<int>(positions.size(), -1);
  }
  std::vector<int> indices;
  for (int i = 0; i < positions.size(); i++) {
    indices.emplace_back(addVertex(positions[i], normals[i]));
  }
  return indices;
}

int LineMesh::connectVertices(int vertex1, int vertex2) {
  return connectVertices(Eigen::Array2i(vertex1, vertex2));
}

int LineMesh::connectVertices(Eigen::Array2i segment) {
  if (!isVertexActive(segment[0]) || !isVertexActive(segment[1]))
    return -1;
  int segId = hasConnection(segment[0], segment[1]);
  if (segId < 0) {
    if (_segmentIdQueue.empty()) {
      segId = _segments.size();
      _segments.emplace_back(0, 0);
      _activeSegment.emplace_back(true);
    } else {
      segId = _segmentIdQueue.front();
      _segmentIdQueue.pop();
    }
    _segments[segId] = segment;
    _vertSegments[segment[0]].emplace_back(segId);
    _vertSegments[segment[1]].emplace_back(segId);
  }
  _activeSegment[segId] = true;
  return segId;
}

std::vector<int>
LineMesh::connectVertices(std::vector<Eigen::Array2i> segments) {
  std::vector<int> indices;
  for (auto segment : segments) {
    int segId = connectVertices(segment);
    if (segId >= 0)
      indices.emplace_back(segId);
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
  if (vertexId < 0 || !_activeVertices[vertexId])
    return 0;
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

void LineMesh::removeVertex(int vertexId) { _activeVertices[vertexId] = false; }

void LineMesh::disconnectVertices(int vertex1, int vertex2) {
  // Find which segment is between those two vertices
  auto segId = hasConnection(vertex1, vertex2);
  if (segId < 0)
    return;

  // In the list of segmetns of each vertex, look for that to be deleted
  std::vector<int> vIds(2);
  vIds[0] = vertex1;
  vIds[1] = vertex2;
  for (auto vId : vIds) {
    auto &segments = getVertexSegments(vId);
    int segPosId = 0;
    while (segments[segPosId] != segId && segPosId < segments.size()) {
      segPosId++;
    }
    if (segPosId < segments.size())
      segments.erase(segments.begin() + segPosId);
  }
  _segmentIdQueue.push(segId);
  _activeSegment[segId] = false;
}

bool LineMesh::isSegmentActive(int segId) { return _activeSegment[segId]; }

bool LineMesh::isVertexActive(int vertexId) {
  return _activeVertices[vertexId];
}

bool LineMesh::isCorner(int vertexId) { return _isCorner[vertexId]; }

Eigen::Vector2d LineMesh::getSegmentNormal(int segId) {
  auto vertices = getSegmentVertices(segId);
  Eigen::Vector2d tangent =
      (getVertexPosition(vertices[1]) - getVertexPosition(vertices[0]));
  return Eigen::Vector2d(tangent[1], -tangent[0]);
}

} // namespace Ramuh