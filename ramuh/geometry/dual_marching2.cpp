#include <geometry/dual_marching.h>
#include <iostream>

namespace Ramuh {

DualMarching2::DualMarching2() { _resolution = 16; }

DualMarching2::DualMarching2(int resolution) { _resolution = resolution; }

Eigen::Array2d
DualMarching2::evaluateSquare(std::pair<int, int> pointIndices,
                              std::vector<Eigen::Array2d> normalLocation,
                              std::vector<Eigen::Vector2d> normals) {

  Eigen::MatrixXd A(normals.size(), 2);
  Eigen::VectorXd b(normals.size());

  for (int i = 0; i < normals.size(); i++) {
    A.row(i) << normals[i][0], normals[i][1];
    b[i] = normals[i].dot(normalLocation[i].matrix());
  }
//   std::cerr << A << std::endl;
//   std::cerr << b << std::endl;
  Eigen::Vector2d x =
      (A.transpose() * A).colPivHouseholderQr().solve(A.transpose() * b);

  _points[pointIndices] = x;
  return x;
}

} // namespace Ramuh