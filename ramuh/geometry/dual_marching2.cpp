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

  //

  return reconstruct(connections);
}

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
  // Create auxiliary points to create slim face, so object can be displayed as
  // a 2D lines
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
  Eigen::Vector2d faceDir = (_points[ids[0]] - _points[ids[1]]).matrix();
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

void DualMarching2::setBaseFolder(std::string folder) { _baseFolder = folder; }

void DualMarching2::resetCounter() { count = 0; }

void DualMarching2::resetCounter(int value) { count = value; }

} // namespace Ramuh