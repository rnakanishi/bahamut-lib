#include <geometry/dual_marching.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace Ramuh {

DualMarching3::DualMarching3() : DualMarching3(Eigen::Array3i(16, 16, 16)) {}

DualMarching3::DualMarching3(Eigen::Array3i resolution) {
  _resolution = resolution;
}

Eigen::Array3i DualMarching3::convertKey(int index) {
  Eigen::Array3i indices(3);
  indices[2] = index / (_resolution[0] * _resolution[1]);
  indices[1] = (index % (_resolution[0] * _resolution[1])) / _resolution[0];
  indices[0] = (index % (_resolution[0] * _resolution[1])) % _resolution[0];
  return indices;
}

int DualMarching3::convertKey(Eigen::Array3i index) {
  int i, j, k;
  i = index[0];
  j = index[1];
  k = index[2];
  return k * _resolution[0] * _resolution[1] + j * _resolution[0] + i;
}

int DualMarching3::convertKey(int i, int j, int k) {
  return k * _resolution[0] * _resolution[1] + j * _resolution[0] + i;
}

Eigen::Array3d
DualMarching3::evaluateCube(Eigen::Array3i pointIndices,
                            std::vector<Eigen::Array3d> normalLocation,
                            std::vector<Eigen::Vector3d> normals) {
  return evaluateCube(pointIndices, normalLocation, normals,
                      BoundingBox3(Eigen::Array3d(-1e8, -1e8, -1e8),
                                   Eigen::Array3d(1e8, 1e8, 1e8)));
}

Eigen::Array3d DualMarching3::evaluateCube(
    Eigen::Array3i cellIndex, std::vector<Eigen::Array3d> normalLocation,
    std::vector<Eigen::Vector3d> normals, BoundingBox3 cubeLimits) {
  if (_idMap.find(convertKey(cellIndex)) != _idMap.end())
    return _points[_idMap[convertKey(cellIndex)]];

  int nsize = normals.size();
  Eigen::MatrixXd A(normals.size(), 3);
  Eigen::VectorXd b(normals.size());
  Eigen::MatrixXd Full(normals.size(), 4);
  Eigen::Vector3d normalAvg(0, 0, 0);
  Eigen::Array3d posAvg(0, 0, 0);

  for (int i = 0; i < normals.size(); i++) {
    A.row(i) << normals[i][0], normals[i][1], normals[i][2];
    b[i] = normals[i].dot(normalLocation[i].matrix());
    // b[i] = 0;
    Full.row(i) << A.row(i), b[i];

    posAvg += normalLocation[i];
    normalAvg += normals[i];
  }
  // posAvg = cubeLimits.center();
  posAvg /= normalLocation.size();

  // Eigen::Vector3d x = (A).colPivHouseholderQr().solve(b);
  // (A.transpose() * A).colPivHouseholderQr().solve(A.transpose() * b);
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(Full.rows(), Full.cols());
  qr.compute(Full);
  Eigen::MatrixXd Q = qr.householderQ();
  Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();

  A = R.block<3, 3>(0, 0);
  b = R.block<3, 1>(0, 3);

  // Pseudo Inverse Computation
  auto svd =
      (A.adjoint() * A).jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
  auto &singularValues = svd.singularValues();
  Eigen::MatrixXd dInv(3, 3);
  dInv.setZero();
  for (size_t i = 0; i < 3; i++) {
    if (singularValues(i) > 0.1)
      dInv(i, i) = 1. / singularValues(i);
    else
      dInv(i, i) = 0.;
  }
  Eigen::Matrix3d pseudoInv = svd.matrixV() * dInv * svd.matrixU().adjoint();

  Eigen::Vector3d x;
  x = pseudoInv * (A.adjoint() * b - A.adjoint() * A * posAvg.matrix());
  x = x + posAvg.matrix();

  if (!cubeLimits.contains(x)) {
    // std::cerr << "Clamp: " << x.transpose() << " -> " << posAvg.transpose()
    // << std::endl;
    x = cubeLimits.clamp(x);
    // x = posAvg;
  }

  int currIndex = _points.size();
  _idMap[convertKey(cellIndex)] = currIndex;
  _points.emplace_back(x);
  _normals.emplace_back(normalAvg);

  return x;
}

bool DualMarching3::_consistentNormals(std::vector<int> ids) {
  Eigen::Vector3d avgNormal;

  for (auto id : ids)
    avgNormal += _normals[id];
  avgNormal /= ids.size();

  Eigen::Vector3d base01 = (_points[ids[1]] - _points[ids[0]]).matrix();
  Eigen::Vector3d base02 = (_points[ids[2]] - _points[ids[0]]).matrix();
  Eigen::Vector3d base03 = (_points[ids[3]] - _points[ids[0]]).matrix();
  Eigen::Vector3d base10 = (_points[ids[0]] - _points[ids[1]]).matrix();
  Eigen::Vector3d base12 = (_points[ids[2]] - _points[ids[1]]).matrix();
  Eigen::Vector3d base13 = (_points[ids[3]] - _points[ids[1]]).matrix();

  Eigen::Vector3d cross012 = base01.cross(base02);
  Eigen::Vector3d cross023 = base02.cross(base03);
  Eigen::Vector3d cross123 = base12.cross(base13);
  Eigen::Vector3d cross130 = base13.cross(base10);

  Eigen::Vector3d avgCross = (cross012 + cross023 + cross123 + cross130) / 4;
  avgCross.normalize();
  avgNormal.normalize();

  if (avgCross.dot(avgNormal) < 0)
    return false;
  return true;
}

void DualMarching3::reconstruct() {

  reconstruct(std::vector<std::pair<int, int>>());
}

void DualMarching3::reconstruct(std::vector<std::pair<int, int>> connections) {
  _buildConnectionMap(connections);

  std::ofstream file;
  static int count = 0;
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

  int id = 1;
  for (int i = 0; i < _points.size(); i++) {
    file << "v " << _points[i][0] << " " << _points[i][1] << " "
         << _points[i][2] << std::endl;
  }
  for (int i = 0; i < _normals.size(); i++) {
    file << "vn " << _normals[i][0] << " " << _normals[i][1] << " "
         << _normals[i][2] << std::endl;
  }

  int fCount = 0;
  // FIXME: When eight surface boundary cells are near, they merge faces
  // wrongly. The outside surface is correct, but inner faces are created too
  for (int i = 0; i < _resolution[0]; i++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int k = 0; k < _resolution[2]; k++) {
        if (_idMap.find(convertKey(i, j, k)) != _idMap.end()) {
          // Verify all directions
          // Right top
          if (_idMap.find(convertKey(i + 1, j, k)) != _idMap.end() &&
              _idMap.find(convertKey(i + 1, j + 1, k)) != _idMap.end() &&
              _idMap.find(convertKey(i, j + 1, k)) != _idMap.end()) {
            if (_hasConnection(Eigen::Array3i(i, j, k),
                               Eigen::Array3i(i + 1, j, k)) &&
                _hasConnection(Eigen::Array3i(i + 1, j, k),
                               Eigen::Array3i(i + 1, j + 1, k)) &&
                _hasConnection(Eigen::Array3i(i + 1, j + 1, k),
                               Eigen::Array3i(i, j + 1, k)) &&
                _hasConnection(Eigen::Array3i(i, j + 1, k),
                               Eigen::Array3i(i, j, k))) {

              fCount++;
              std::vector<int> pIds;
              pIds.emplace_back(_idMap[convertKey(i, j, k)]);
              pIds.emplace_back(_idMap[convertKey(i + 1, j, k)]);
              pIds.emplace_back(_idMap[convertKey(i + 1, j + 1, k)]);
              pIds.emplace_back(_idMap[convertKey(i, j + 1, k)]);
              if (!_consistentNormals(pIds))
                std::swap(pIds[0], pIds[2]);
              file << "f ";
              for (auto id : pIds) {
                file << id + 1 << " ";
              }
              file << std::endl;
            }
          }
          // Front right
          if (_idMap.find(convertKey(i, j, k + 1)) != _idMap.end() &&
              _idMap.find(convertKey(i + 1, j, k + 1)) != _idMap.end() &&
              _idMap.find(convertKey(i + 1, j, k)) != _idMap.end()) {
            if (_hasConnection(Eigen::Array3i(i, j, k),
                               Eigen::Array3i(i, j, k + 1)) &&
                _hasConnection(Eigen::Array3i(i, j, k + 1),
                               Eigen::Array3i(i + 1, j, k + 1)) &&
                _hasConnection(Eigen::Array3i(i + 1, j, k + 1),
                               Eigen::Array3i(i + 1, j, k)) &&
                _hasConnection(Eigen::Array3i(i + 1, j, k),
                               Eigen::Array3i(i, j, k))) {
              fCount++;
              std::vector<int> pIds;
              pIds.emplace_back(_idMap[convertKey(i, j, k)]);
              pIds.emplace_back(_idMap[convertKey(i, j, k + 1)]);
              pIds.emplace_back(_idMap[convertKey(i + 1, j, k + 1)]);
              pIds.emplace_back(_idMap[convertKey(i + 1, j, k)]);
              if (!_consistentNormals(pIds))
                std::swap(pIds[0], pIds[2]);
              file << "f ";
              for (auto id : pIds) {
                file << id + 1 << " ";
              }
              file << std::endl;
            }
          }
          // top front
          if (_idMap.find(convertKey(i, j + 1, k)) != _idMap.end() &&
              _idMap.find(convertKey(i, j + 1, k + 1)) != _idMap.end() &&
              _idMap.find(convertKey(i, j, k + 1)) != _idMap.end()) {
            if (_hasConnection(Eigen::Array3i(i, j, k),
                               Eigen::Array3i(i, j + 1, k)) &&
                _hasConnection(Eigen::Array3i(i, j + 1, k),
                               Eigen::Array3i(i, j + 1, k + 1)) &&
                _hasConnection(Eigen::Array3i(i, j + 1, k + 1),
                               Eigen::Array3i(i, j, k + 1)) &&
                _hasConnection(Eigen::Array3i(i, j, k + 1),
                               Eigen::Array3i(i, j, k))) {
              fCount++;
              std::vector<int> pIds;
              pIds.emplace_back(_idMap[convertKey(i, j, k)]);
              pIds.emplace_back(_idMap[convertKey(i, j + 1, k)]);
              pIds.emplace_back(_idMap[convertKey(i, j + 1, k + 1)]);
              pIds.emplace_back(_idMap[convertKey(i, j, k + 1)]);
              if (!_consistentNormals(pIds))
                std::swap(pIds[0], pIds[2]);
              file << "f ";
              for (auto id : pIds) {
                file << id + 1 << " ";
              }
              file << std::endl;
            }
          }
        }
      }
  std::cerr << "File written: " << filename.str();
  std::cerr << ": " << _points.size() << " vertices, " << fCount << " faces.\n";
} // namespace Ramuh

std::vector<Eigen::Array3d> &DualMarching3::getPoints() { return _points; }

std::vector<Eigen::Vector3d> &DualMarching3::getNormals() { return _normals; }

std::map<int, int> &DualMarching3::getIdMap() { return _idMap; }

void DualMarching3::merge(DualMarching3 cube) {
  auto &points = cube.getPoints();
  auto &normals = cube.getNormals();
  auto &map = cube.getIdMap();
  int nverts = _points.size();

  _points.insert(_points.end(), points.begin(), points.end());
  _normals.insert(_normals.end(), normals.begin(), normals.end());

  for (auto &&point : map) {
    point.second += nverts;
  }
  _idMap.insert(map.begin(), map.end());
}

bool DualMarching3::_hasConnection(Eigen::Array3i tuple1,
                                   Eigen::Array3i tuple2) {
  if (_connections.empty())
    return true;

  int vertex1 = _idMap[convertKey(tuple1)];
  for (auto connection : _connections[vertex1]) {
    if (connection == _idMap[convertKey(tuple2)])
      return true;
  }
  return false;
}

void DualMarching3::_buildConnectionMap(
    std::vector<std::pair<int, int>> connections) {

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

void DualMarching3::setBaseFolder(std::string folder) { _baseFolder = folder; }

} // namespace Ramuh