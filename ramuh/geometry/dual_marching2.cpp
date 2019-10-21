#include <geometry/dual_marching.h>
#include <iostream>

namespace Ramuh {

DualMarching2::DualMarching2() { _resolution = 16; }

DualMarching2::DualMarching2(Eigen::Array2i resolution) {
  _resolution = resolution;
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

  int nsize = normals.size();
  Eigen::MatrixXd A(normals.size() + 2, 2);
  Eigen::VectorXd b(normals.size() + 2);
  Eigen::Vector2d normalAvg(0, 0);
  Eigen::Array2d posAvg(0, 0);

  for (int i = 0; i < nsize; i++) {
    A.row(i) << normals[i][0], normals[i][1];
    b[i] = normals[i].dot(normalLocation[i].matrix());
    b[i] = 0;
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
  b[nsize + 0] = Eigen::Vector2d(1e-1, 0.).dot(posAvg.matrix());
  b[nsize + 1] = Eigen::Vector2d(0., 1e-1).dot(posAvg.matrix());

  // Check boundaries so the point remains inside, even if the ouside result is
  // correct

  Eigen::Vector2d x = (A).colPivHouseholderQr().solve(b);
  // (A.transpose() * A).colPivHouseholderQr().solve(A.transpose() * b);

  if (!squareLimits.contains(x)) {
    // x = cubeLimits.clamp(x);
    x = posAvg;
    // TODO: implement constrained QEF solver
  }

  int currIndex = _points.size();
  _idMap[convertKey(cellIndex)] = currIndex;
  _points.emplace_back(x);
  _normals.emplace_back(normalAvg);

  return x;
}

} // namespace Ramuh