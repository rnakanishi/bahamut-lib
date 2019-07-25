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
  _points[nPoints] = x;
  _normals[nPoints] = normalAvg.normalized();

  return x;
}

void DualMarching3::reconstruct() {
  std::ofstream file;
  static int count = 0;
  std::stringstream filename;
  filename << "results/marching/dualcubes_" << std::setfill('0') << std::setw(4)
           << count++ << ".obj";
  file.open(filename.str().c_str(), std::ofstream::out);

  int id = 1;
  for (std::map<int, Eigen::Array3d>::iterator it = _points.begin();
       it != _points.end(); it++) {
    file << "v " << it->second[0] << " " << it->second[1] << " "
         << it->second[2] << std::endl;
  }
  for (std::map<int, Eigen::Vector3d>::iterator it = _normals.begin();
       it != _normals.end(); it++) {
    file << "vn " << it->second[0] << " " << it->second[1] << " "
         << it->second[2] << std::endl;
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
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i + 1, j, k)] << " "
                 << _idMap[std::make_tuple(i + 1, j + 1, k)] << " "
                 << _idMap[std::make_tuple(i, j + 1, k)] << " " << std::endl;
          }
          // Bottom right
          if (_idMap.find(std::make_tuple(i, j - 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j - 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j, k)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j - 1, k)] << " "
                 << _idMap[std::make_tuple(i + 1, j - 1, k)] << " "
                 << _idMap[std::make_tuple(i + 1, j, k)] << " " << std::endl;
          }
          // Left bottom
          if (_idMap.find(std::make_tuple(i - 1, j, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j - 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j - 1, k)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i - 1, j, k)] << " "
                 << _idMap[std::make_tuple(i - 1, j - 1, k)] << " "
                 << _idMap[std::make_tuple(i, j - 1, k)] << " " << std::endl;
          }
          // Top left
          if (_idMap.find(std::make_tuple(i, j + 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j + 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j, k)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j + 1, k)] << " "
                 << _idMap[std::make_tuple(i - 1, j + 1, k)] << " "
                 << _idMap[std::make_tuple(i - 1, j, k)] << " " << std::endl;
          }
          // Right back
          if (_idMap.find(std::make_tuple(i + 1, j, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j, k - 1)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i + 1, j, k)] << " "
                 << _idMap[std::make_tuple(i + 1, j, k - 1)] << " "
                 << _idMap[std::make_tuple(i, j, k - 1)] << " " << std::endl;
          }
          // Back left
          if (_idMap.find(std::make_tuple(i, j, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i + 1, j, k)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j, k - 1)] << " "
                 << _idMap[std::make_tuple(i + 1, j, k - 1)] << " "
                 << _idMap[std::make_tuple(i + 1, j, k)] << " " << std::endl;
          }
          // Left front
          if (_idMap.find(std::make_tuple(i - 1, j, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i - 1, j, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j, k + 1)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i - 1, j, k)] << " "
                 << _idMap[std::make_tuple(i - 1, j, k + 1)] << " "
                 << _idMap[std::make_tuple(i, j, k + 1)] << " " << std::endl;
          }
          // Front right
          if (_idMap.find(std::make_tuple(i, j, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j, k + 1)] << " "
                 << _idMap[std::make_tuple(i, j + 1, k + 1)] << " "
                 << _idMap[std::make_tuple(i, j + 1, k)] << " " << std::endl;
          }
          // back top
          if (_idMap.find(std::make_tuple(i, j, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j, k - 1)] << " "
                 << _idMap[std::make_tuple(i, j + 1, k - 1)] << " "
                 << _idMap[std::make_tuple(i, j + 1, k)] << " " << std::endl;
          }
          // top front
          if (_idMap.find(std::make_tuple(i, j + 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j + 1, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j, k + 1)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j, k)] << " " << std::endl;
          }
          // Front bottom
          if (_idMap.find(std::make_tuple(i, j, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j - 1, k + 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j - 1, k)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j, k + 1)] << " "
                 << _idMap[std::make_tuple(i, j - 1, k + 1)] << " "
                 << _idMap[std::make_tuple(i, j - 1, k)] << " " << std::endl;
          }
          // Bottom back
          if (_idMap.find(std::make_tuple(i, j - 1, k)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j - 1, k - 1)) != _idMap.end() &&
              _idMap.find(std::make_tuple(i, j, k - 1)) != _idMap.end()) {
            file << "f " << _idMap[std::make_tuple(i, j, k)] << " "
                 << _idMap[std::make_tuple(i, j - 1, k)] << " "
                 << _idMap[std::make_tuple(i, j - 1, k - 1)] << " "
                 << _idMap[std::make_tuple(i, j, k - 1)] << " " << std::endl;
          }
        }
      }
  std::cerr << "File written: " << filename.str() << std::endl;
}

} // namespace Ramuh