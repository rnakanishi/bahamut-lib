#include <geometry/dual_marching.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace Ramuh {

DualMarching2::DualMarching2() : DualMarching2(Eigen::Array2i(16, 16)) {}

DualMarching2::DualMarching2(Eigen::Array2i resolution) {
  _resolution = resolution;
  resetCounter();
  clear();
}

void DualMarching2::clear() {
  _idMap.clear();
  _points.clear();
  _normals.clear();
  _connections.clear();
}

int DualMarching2::convertKey(int i, int j) { return j * _resolution[0] + i; }

int DualMarching2::convertKey(Eigen::Array2i index) {
  return convertKey(index[0], index[1]);
}

Eigen::Array2i DualMarching2::convertKey(int id) {
  Eigen::Array2i index;
  index[0] = id % (_resolution[0]);
  index[1] = id / _resolution[0];
  return index;
}

std::vector<Eigen::Array2d> &DualMarching2::getPoints() { return _points; }

std::vector<Eigen::Vector2d> &DualMarching2::getNormals() { return _normals; }

std::map<int, int> &DualMarching2::getIdMap() { return _idMap; }

Eigen::Array2d
DualMarching2::evaluateSquare(Eigen::Array2i pointIndices,
                              std::vector<Eigen::Array2d> normalLocation,
                              std::vector<Eigen::Vector2d> normals) {
  return evaluateSquare(
      pointIndices, normalLocation, normals,
      BoundingBox2(Eigen::Array2d(-1e8, -1e8), Eigen::Array2d(1e8, 1e8)));
}

Eigen::Array2d DualMarching2::evaluateSquare(
    Eigen::Array2i cellIndex, std::vector<Eigen::Array2d> normalLocation,
    std::vector<Eigen::Vector2d> normals, BoundingBox2 squareLimits) {
  if (_idMap.find(convertKey(cellIndex)) != _idMap.end())
    return _points[_idMap[convertKey(cellIndex)]];

  int validNormals = 0;
  int nsize = normals.size();
  Eigen::MatrixXd A(normals.size(), 2);
  Eigen::VectorXd b(normals.size());
  Eigen::MatrixXd Full(normals.size(), 3);
  Eigen::Vector2d normalAvg(0, 0);
  Eigen::Array2d posAvg(0, 0);
  Eigen::Vector2d x;

  if (nsize > 1) {
    // If only one point is given, minimization doesn`t have to be computed
    for (int i = 0; i < nsize; i++) {
      A.row(i) << normals[i][0], normals[i][1];
      b[i] = normals[i].dot(normalLocation[i].matrix());
      Full.row(i) << A.row(i), b[i];

      if (normalLocation[i].matrix().norm() > 1e-6) {
        validNormals++;
        posAvg += normalLocation[i];
        normalAvg += normals[i];
      }
    }
    // posAvg = cubeLimits.center();
    posAvg /= validNormals;
    // TODO: check if normals make a corner
    normalAvg /= validNormals;

    // QR minimization
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(Full.rows(), Full.cols());
    qr.compute(Full);
    Eigen::MatrixXd Q = qr.householderQ();
    Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
    A = R.block<2, 2>(0, 0);
    b = R.block<2, 1>(0, 2);

    // Pseudo Inverse Computation
    auto svd =
        (A.adjoint() * A).jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto &singularValues = svd.singularValues();
    Eigen::MatrixXd dInv(2, 2);
    dInv.setZero();
    for (size_t i = 0; i < 2; i++) {
      if (singularValues(i) > 0.1)
        dInv(i, i) = 1. / singularValues(i);
      else
        dInv(i, i) = 0.;
    }
    Eigen::Matrix2d pseudoInv = svd.matrixV() * dInv * svd.matrixU().adjoint();

    x = pseudoInv * (A.adjoint() * b - A.adjoint() * A * posAvg.matrix());
    x = x + posAvg.matrix();

    // Check boundaries so the point remains inside, even if the ouside result
    // is correct
    if (!squareLimits.contains(x)) {
      // x = cubeLimits.clamp(x);
      x = posAvg;
    }
  } else {
    x = normalLocation[0];
    normalAvg = normals[0];
  }
  int currIndex = _points.size();
  _idMap[convertKey(cellIndex)] = currIndex;
  _points.emplace_back(x);
  _normals.emplace_back(normalAvg);

  return x;
}

Ramuh::LineMesh DualMarching2::reconstruct() {
  // Build connections
  std::vector<std::pair<int, int>> connections;
  std::vector<bool> visited(_points.size(), false);
  Ramuh::LineMesh mesh;
  std::queue<int> simpleConn, tripleConn;

  // Add vertices and their connections
  _createSimpleConnections(mesh);

  // Look for spurious vertices and fix them
  for (auto cell : _idMap) {
    int vertexId = cell.second;
    if (!mesh.isVertexActive(vertexId))
      continue;
    auto ij = convertKey(cell.first);

    // If vertex connection is different from 2, then it is spurious
    int nConnections = mesh.getNumberOfConnections(vertexId);
    if (nConnections < 2)
      simpleConn.push(cell.first);
    else if (nConnections > 2)
      tripleConn.push(cell.first);
  }

  // Process all single connections first
  while (!simpleConn.empty()) {
    int cell = simpleConn.front();
    simpleConn.pop();
    int vertexId = _idMap[cell];
    if (mesh.getNumberOfConnections(vertexId) == 2)
      continue;

    auto neighbors = mesh.getAdjacentVertices(vertexId);
    int neighborId = -1;
    // First try to find the other vertex that is also missing a connection
    bool connected = false;
    auto ij = convertKey(cell);
    if (neighbors.size() > 0)
      neighborId = neighbors[0];
    if (mesh.getNumberOfConnections(neighborId) < 3) {
      std::vector<int> toConnect;
      for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
          if (i == 0 && j == 0)
            continue;
          int neighCell = convertKey(ij[0] + i, ij[1] + j);
          if (_idMap.find(neighCell) != _idMap.end()) {
            if (mesh.getNumberOfConnections(_idMap[neighCell]) <= 1) {
              connected = true;
              // Simple case: both need to be connected
              if (checkOrientation(mesh, vertexId, _idMap[neighCell]))
                mesh.connectVertices(vertexId, _idMap[neighCell]);
              else
                mesh.connectVertices(_idMap[neighCell], vertexId);
            } else if (mesh.getNumberOfConnections(vertexId) <= 1 &&
                       mesh.getNumberOfConnections(_idMap[neighCell]) == 2) {
              // Isolation: Connect both and the neighbor vertex will be
              // analysed
              toConnect.emplace_back(neighCell);
            }
          }
        }
      }
      if (!connected) {
        for (auto cellNeigh : toConnect) {
          int neighId = _idMap[cellNeigh];
          if (mesh.hasConnection(vertexId, neighId) < 0) {
            tripleConn.push(cellNeigh);
            // if (checkNormalDirection(vertexId, neighId))
            if (checkOrientation(mesh, vertexId, neighId))
              mesh.connectVertices(vertexId, neighId);
            else
              mesh.connectVertices(neighId, vertexId);
          }
        }
        if (mesh.getNumberOfConnections(vertexId) == 3)
          tripleConn.push(cell);
      }
    } else {
      // Not able to find any missing neighbor vertex so niehgbor is probably a
      // hub. In this case, its connections are given to the current vertex and
      // the neighbor is disconnected
      auto neighbor = mesh.getAdjacentVertices(vertexId);
      int neighborConnections = mesh.getNumberOfConnections(neighbor[0]);
      auto neighSegments = mesh.getVertexSegments(neighbor[0]);
      auto segmentId = mesh.hasConnection(neighbor[0], vertexId);
      for (auto segment : neighSegments) {
        if (segment != segmentId) {
          auto vertices = mesh.getSegmentVertices(segment);
          if (vertices[0] != neighbor[0]) {
            // if (checkNormalDirection(vertices[0], vertexId))
            if (checkOrientation(mesh, vertices[0], vertexId))
              mesh.connectVertices(vertices[0], vertexId);
            else
              mesh.connectVertices(vertexId, vertices[0]);
            mesh.disconnectVertices(neighbor[0], vertices[0]);
          } else {
            // if (checkNormalDirection(vertices[0], vertexId))
            if (checkOrientation(mesh, vertices[1], vertexId))
              mesh.connectVertices(vertices[1], vertexId);
            else
              mesh.connectVertices(vertexId, vertices[1]);
            mesh.disconnectVertices(neighbor[0], vertices[1]);
          }
        }
      }
      mesh.disconnectVertices(neighbor[0], vertexId);
      mesh.removeVertex(neighbor[0]);
    }
  }
  while (!tripleConn.empty()) {
    int cell = tripleConn.front();
    tripleConn.pop();
    int vertexId = _idMap[cell];
    // Check if stiil has triple connections
    if (mesh.getNumberOfConnections(vertexId) <= 2)
      continue;

    // Check the neighbors
    int lonelyVertex = -1;
    std::vector<int> tripleConnectionIds;
    auto neighbors = mesh.getAdjacentVertices(vertexId);
    for (auto neighbor : neighbors) {
      int nConn = mesh.getNumberOfConnections(neighbor);
      if (nConn == 1) {
        lonelyVertex = neighbor;
      } else if (nConn == 3) {
        tripleConnectionIds.emplace_back(neighbor);
      }
    }
    // Found on neighbor with only one connection: ignore. This case will be
    // treated by that lonely vertex
    if (lonelyVertex >= 0 && tripleConnectionIds.empty()) {
      auto pos1 = mesh.getVertexPosition(vertexId);
      auto pos2 = mesh.getVertexPosition(lonelyVertex);
      if ((pos1 - pos2).matrix().norm() < 1e-5) {
        mesh.disconnectVertices(lonelyVertex, vertexId);
        mesh.removeVertex(lonelyVertex);
      }
      continue;
    }

    // If only one neighbor with three connections is found, than remove the
    // connection between them
    // If more then one neighbor with more than two connections is found,
    // then the common edge has to be found. Disconnect that edge and remove
    // the problematic vertex
    if (tripleConnectionIds.size() == 1) {
      auto neighbor = tripleConnectionIds[0];
      mesh.disconnectVertices(vertexId, neighbor);
    } else if (tripleConnectionIds.size() == 2) {
      if (mesh.hasConnection(tripleConnectionIds[0], tripleConnectionIds[1]))
        mesh.disconnectVertices(tripleConnectionIds[0], tripleConnectionIds[1]);
      else {
      }
    } else if (tripleConnectionIds.size() == 3) {
    }
  }

  // _writeMesh(mesh);
  _mesh = mesh;
  return mesh;
} // namespace Ramuh

Ramuh::LineMesh
DualMarching2::reconstruct(std::vector<std::pair<int, int>> connections) {
  _buildConnectionMap(connections);

  std::ofstream file;
  std::stringstream filename;
  filename << _baseFolder << std::setfill('0') << std::setw(4) << count++
           << ".obj";
  auto fullpath = filename.str();
  file.open(fullpath.c_str(), std::ofstream::out);
  if (!file.is_open()) {
    std::cerr << "\033[1;31m[ X ]\033[0m Failed opening file " << filename.str()
              << std::endl;
    return Ramuh::LineMesh();
  }

  Ramuh::LineMesh mesh;
  int id = 1;
  int realPointCount = _points.size();
  for (int i = 0; i < _points.size(); i++) {
    file << "v " << _points[i][0] << " " << _points[i][1] << " 0" << std::endl;
    mesh.addVertex(_points[i]);
  }
  // Create auxiliary points to create slim face, so object can be displayed
  // as a 2D lines
  for (int i = 0; i < _points.size(); i++) {
    file << "v " << _points[i][0] << " " << _points[i][1] << " 0.001"
         << std::endl;
  }

  int fCount = 0;
  // FIXME: When eight surface boundary cells are near, they merge faces
  // wrongly. The outside surface is correct, but inner faces are created too
  for (int i = 0; i < _resolution[0]; i++)
    for (int j = 0; j < _resolution[1]; j++)
      if (_idMap.find(convertKey(i, j)) != _idMap.end()) {
        // Verify all directions
        // Right
        if (_idMap.find(convertKey(i + 1, j)) != _idMap.end()) {
          if (_hasConnection(Eigen::Array2i(i, j), Eigen::Array2i(i + 1, j))) {
            fCount++;
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[convertKey(i, j)]);
            pIds.emplace_back(_idMap[convertKey(i + 1, j)]);
            pIds.emplace_back(_idMap[convertKey(i + 1, j)] + realPointCount);
            pIds.emplace_back(_idMap[convertKey(i, j)] + realPointCount);
            file << "f ";
            if (!_consistentNormals(pIds)) {
              std::swap(pIds[0], pIds[1]);
              std::swap(pIds[2], pIds[3]);
            }
            for (auto id : pIds) {
              file << id + 1 << " ";
            }
            mesh.connectVertices(Eigen::Array2i(pIds[0], pIds[1]));
            file << std::endl;
          }
        }
        // Top
        if (_idMap.find(convertKey(i, j + 1)) != _idMap.end()) {
          if (_hasConnection(Eigen::Array2i(i, j), Eigen::Array2i(i, j + 1))) {
            fCount++;
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[convertKey(i, j)]);
            pIds.emplace_back(_idMap[convertKey(i, j + 1)]);
            pIds.emplace_back(_idMap[convertKey(i, j + 1)] + realPointCount);
            pIds.emplace_back(_idMap[convertKey(i, j)] + realPointCount);
            file << "f ";
            if (!_consistentNormals(pIds)) {
              std::swap(pIds[0], pIds[1]);
              std::swap(pIds[2], pIds[3]);
            }
            for (auto id : pIds) {
              file << id + 1 << " ";
            }
            mesh.connectVertices(Eigen::Array2i(pIds[0], pIds[1]));
            file << std::endl;
          }
        }
      }

  std::cerr << "File written: " << filename.str();
  std::cerr << ": " << _points.size() << " vertices, " << fCount << " faces.\n";
  return mesh;
}

bool DualMarching2::_consistentNormals(std::vector<int> ids) {
  Eigen::Vector2d faceDir = (_points[ids[1]] - _points[ids[0]]).matrix();
  Eigen::Vector2d faceNormal(faceDir[1], -faceDir[0]);
  if (faceNormal.dot(_normals[ids[0]]) < 0 ||
      faceNormal.dot(_normals[ids[1]]) < 0)
    return false;
  return true;
}

bool DualMarching2::_hasConnection(Eigen::Array2i tuple1,
                                   Eigen::Array2i tuple2) {
  if (_connections.empty())
    return true;

  int vertex1 = _idMap[convertKey(tuple1)];
  for (auto connection : _connections[vertex1]) {
    if (connection == _idMap[convertKey(tuple2)])
      return true;
  }
  return false;
}

void DualMarching2::_buildConnectionMap(
    std::vector<std::pair<int, int>> connections) {

  _connections.clear();
  for (auto connection : connections) {
    int id1 = _idMap[connection.first];
    int id2 = _idMap[connection.second];

    if (_connections.find(id1) == _connections.end()) {
      _connections[id1] = std::vector<int>();
    }
    if (_connections.find(id2) == _connections.end()) {
      _connections[id2] = std::vector<int>();
    }
    _connections[id1].emplace_back(id2);
    _connections[id2].emplace_back(id1);
  }
}

void DualMarching2::merge(DualMarching2 square) {
  auto &points = square.getPoints();
  auto &normals = square.getNormals();
  auto &map = square.getIdMap();
  int nverts = _points.size();

  _points.insert(_points.end(), points.begin(), points.end());
  _normals.insert(_normals.end(), normals.begin(), normals.end());

  for (auto &&point : map) {
    point.second += nverts;
  }
  _idMap.insert(map.begin(), map.end());
}

bool DualMarching2::checkOrientation(LineMesh &mesh, int vertex, int target) {
  Eigen::Array2d origin, ending;
  origin = mesh.getVertexPosition(vertex);
  ending = mesh.getVertexPosition(target);
  Eigen::Vector2d vector = (ending - origin).matrix().normalized();
  Eigen::Vector2d normal(vector[1], -vector[0]);
  Eigen::Vector2d vNormal = _normals[vertex];
  if (vNormal.dot(normal) > 0)
    return true;
  return false;
}

bool DualMarching2::checkNormalDirection(int vert1Id, int vert2Id) {
  Eigen::Vector2d normal1, normal2;
  normal1 = _normals[vert1Id];
  normal2 = _normals[vert2Id];

  if (normal1.dot(normal2) < 0)
    return false;
  return true;
}

void DualMarching2::_createSimpleConnections(LineMesh &mesh) {

  mesh.addVertices(_points, _normals);

  // for all cells, evaluate 4-neighborhood and connect with all cells that
  // have a vertex
  for (auto cell : _idMap) {
    int cellId = cell.first;
    auto ij = convertKey(cellId);
    std::vector<int> neighsToConnect;

    // Check 4-neighborhood
    if (ij[0] > 0 &&
        _idMap.find(convertKey(ij[0] - 1, ij[1])) != _idMap.end()) {
      neighsToConnect.emplace_back(convertKey(ij[0] - 1, ij[1]));
    }
    if (ij[0] < _resolution[0] - 1 &&
        _idMap.find(convertKey(ij[0] + 1, ij[1])) != _idMap.end()) {
      neighsToConnect.emplace_back(convertKey(ij[0] + 1, ij[1]));
    }
    if (ij[1] > 0 &&
        _idMap.find(convertKey(ij[0], ij[1] - 1)) != _idMap.end()) {
      neighsToConnect.emplace_back(convertKey(ij[0], ij[1] - 1));
    }
    if (ij[1] < _resolution[1] - 1 &&
        _idMap.find(convertKey(ij[0], ij[1] + 1)) != _idMap.end()) {
      neighsToConnect.emplace_back(convertKey(ij[0], ij[1] + 1));
    }
    for (auto neighId : neighsToConnect) {
      // check if the vertices are not pointing against each other
      if (checkNormalDirection(_idMap[cellId], _idMap[neighId]))
        // Check face normal orientation
        if (checkOrientation(mesh, _idMap[cellId], _idMap[neighId]))
          mesh.connectVertices(_idMap[cellId], _idMap[neighId]);
        else
          mesh.connectVertices(_idMap[neighId], _idMap[cellId]);
    }
  }
}

void DualMarching2::setBaseFolder(std::string folder) { _baseFolder = folder; }

void DualMarching2::resetCounter() { count = 0; }

void DualMarching2::resetCounter(int value) { count = value; }

void DualMarching2::_writeMesh(LineMesh &mesh) {
  std::ofstream file;
  std::stringstream filename;
  filename << _baseFolder << std::setfill('0') << std::setw(4) << count++
           << ".obj";
  auto fullpath = filename.str();
  file.open(fullpath.c_str(), std::ofstream::out);
  if (!file.is_open()) {
    std::cerr << "\033[1;31m[ X ]\033[0m Failed opening file " << filename.str()
              << std::endl;
    return;
  }

  auto &vertices = mesh.getVerticesList();
  auto &segments = mesh.getSegmentsList();
  int nVertices = vertices.size();
  std::vector<int> inactiveVertices;
  inactiveVertices.clear();
  for (size_t vId = 0; vId < vertices.size(); vId++) {
    auto vertex = vertices[vId];
    if (mesh.isVertexActive(vId))
      file << "v " << vertex[0] << " " << vertex[1] << " 0" << std::endl;
    else {
      inactiveVertices.emplace_back(vId);
      nVertices--;
    }
  }
  for (size_t vId = 0; vId < vertices.size(); vId++) {
    auto vertex = vertices[vId];
    if (mesh.isVertexActive(vId))
      file << "v " << vertex[0] << " " << vertex[1] << " 1" << std::endl;
  }

  for (size_t sId = 0; sId < segments.size(); sId++) {
    if (!mesh.isSegmentActive(sId))
      continue;
    auto segment = segments[sId];
    Eigen::Array2i reduction(0, 0);
    for (auto vertex : inactiveVertices) {
      if (segment[0] >= vertex)
        reduction[0]++;
      if (segment[1] >= vertex)
        reduction[1]++;
    }
    segment = segment - reduction;

    file << "f " << segment[0] + 1 << " " << segment[1] + 1 << " "
         << segment[1] + 1 + nVertices << " " << segment[0] + 1 + nVertices
         << " " << std::endl;
  }
  std::cerr << "File written: " << filename.str();
  std::cerr << ": " << vertices.size() << " vertices, " << segments.size()
            << " faces.\n";
}

} // namespace Ramuh