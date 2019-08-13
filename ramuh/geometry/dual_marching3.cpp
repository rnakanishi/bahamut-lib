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

Eigen::Array3d
DualMarching3::evaluateCube(std::tuple<int, int, int> pointIndices,
                            std::vector<Eigen::Array3d> normalLocation,
                            std::vector<Eigen::Vector3d> normals) {
  return evaluateCube(pointIndices, normalLocation, normals,
                      BoundingBox3(Eigen::Array3d(-1e8, -1e8, -1e8),
                                   Eigen::Array3d(1e8, 1e8, 1e8)));
}
Eigen::Array3d
DualMarching3::evaluateCube(std::tuple<int, int, int> pointIndices,
                            std::vector<Eigen::Array3d> normalLocation,
                            std::vector<Eigen::Vector3d> normals,
                            BoundingBox3 cubeLimits) {
  int nsize = normals.size();
  Eigen::MatrixXd A(normals.size() + 3, 3);
  Eigen::VectorXd b(normals.size() + 3);
  Eigen::Vector3d normalAvg(0, 0, 0);
  Eigen::Array3d posAvg(0, 0, 0);

  for (int i = 0; i < normals.size(); i++) {
    A.row(i) << normals[i][0], normals[i][1], normals[i][2];
    b[i] = normals[i].dot(normalLocation[i].matrix());
    posAvg += normalLocation[i];
    normalAvg += normals[i];
  }
  // posAvg = cubeLimits.center();
  posAvg /= normalLocation.size();
  normalAvg /= normals.size();

  // Bias
  // Adding more vectors so it enforces the new point to be inside the cube
  A.row(nsize + 0) << 1e-1, 0., 0.;
  A.row(nsize + 1) << 0., 1e-1, 0.;
  A.row(nsize + 2) << 0., 0., 1e-1;
  b[nsize + 0] = Eigen::Vector3d(1e-1, 0., 0.).dot(posAvg.matrix());
  b[nsize + 1] = Eigen::Vector3d(0., 1e-1, 0.).dot(posAvg.matrix());
  b[nsize + 2] = Eigen::Vector3d(0., 0., 1e-1).dot(posAvg.matrix());

  // Check boundaries so the point remains inside, even if the ouside result is
  // correct

  Eigen::Vector3d x =
      // (A.transpose() * A).colPivHouseholderQr().solve(A.transpose() * b);
      (A).colPivHouseholderQr().solve(b);

  if (!cubeLimits.contains(x)) {
    // std::cerr << "Limits: " << cubeLimits.min().transpose() << " x "
    //           << cubeLimits.max().transpose();
    // std::cerr << "   point: " << x.transpose() << std::endl;
    x = cubeLimits.clamp(x);
    // TODO: implement constrained QEF solver
  }

  int nPoints = _idMap.size() + 1;
  _idMap[pointIndices] = nPoints;
  _points.emplace_back(x);
  _normals.emplace_back(normalAvg.normalized());

  return x;
}

bool DualMarching3::_consistentNormals(std::vector<int> ids) {
  Eigen::Vector3d base1, base2;

  base1 = (_points[ids[1]] - _points[ids[0]]).matrix();
  base2 = (_points[ids[2]] - _points[ids[1]]).matrix();
  if (base1.cross(base2).dot(_normals[ids[1]]) < 0)
    return false;
  return true;
}

void DualMarching3::reconstruct() {
  std::ofstream file;
  static int count = 0;
  std::stringstream filename;
  filename << "results/marching/dualcubes_" << std::setfill('0') << std::setw(4)
           << count++ << ".obj";
  try {
    file.open(filename.str().c_str(), std::ofstream::out);
  } catch (std::exception e) {
    std::cerr << "Failed opening file " << filename.str() << std::endl;
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

  for (int i = 0; i < _resolution[0]; i++)
    for (int j = 0; j < _resolution[1]; j++)
      for (int k = 0; k < _resolution[2]; k++) {
        if (_idMap.find(std::make_tuple(i, j, k)) != _idMap.end()) {
          // Verify all directions
          // Right top
          if (_idMap.find(std::make_tuple(i + 1, j, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j + 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i + 1, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i + 1, j + 1, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j + 1, k)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Bottom right
          if (_idMap.find(std::make_tuple(i, j - 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j - 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j, k)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j - 1, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i + 1, j - 1, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i + 1, j, k)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Left bottom
          if (_idMap.find(std::make_tuple(i - 1, j, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j - 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j - 1, k)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i - 1, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i - 1, j - 1, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j - 1, k)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Top left
          if (_idMap.find(std::make_tuple(i, j + 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j + 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j, k)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j + 1, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i - 1, j + 1, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i - 1, j, k)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Right back
          if (_idMap.find(std::make_tuple(i + 1, j, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j, k - 1)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i + 1, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i + 1, j, k - 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k - 1)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Back left
          if (_idMap.find(std::make_tuple(i, j, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j, k)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k - 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i - 1, j, k - 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i - 1, j, k)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Left front
          if (_idMap.find(std::make_tuple(i - 1, j, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j, k + 1)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i - 1, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i - 1, j, k + 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k + 1)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Front right
          if (_idMap.find(std::make_tuple(i, j, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j, k)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k + 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i + 1, j, k + 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i + 1, j, k)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // back top
          if (_idMap.find(std::make_tuple(i, j, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k - 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j + 1, k - 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j + 1, k)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // top front
          if (_idMap.find(std::make_tuple(i, j + 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j, k + 1)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j + 1, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j + 1, k + 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k + 1)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Front bottom
          if (_idMap.find(std::make_tuple(i, j, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j - 1, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j - 1, k)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k + 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j - 1, k + 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j - 1, k)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
          // Bottom back
          if (_idMap.find(std::make_tuple(i, j - 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j - 1, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j, k - 1)) != _idMap.end()) {
            std::vector<int> pIds;
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j - 1, k)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j - 1, k - 1)]);
            pIds.emplace_back(_idMap[std::make_tuple(i, j, k - 1)]);
            if (!_consistentNormals(pIds))
              std::swap(pIds[0], pIds[2]);
            file << "f ";
            for (auto id : pIds) {
              file << id << " ";
            }
            file << std::endl;
          }
        }
      }
  std::cerr << "File written: " << filename.str() << std::endl;
} // namespace Ramuh

std::vector<Eigen::Array3d> &DualMarching3::getPoints() { return _points; }

std::vector<Eigen::Vector3d> &DualMarching3::getNormals() { return _normals; }

std::map<std::tuple<int, int, int>, int> &DualMarching3::getIdMap() {
  return _idMap;
}

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

} // namespace Ramuh