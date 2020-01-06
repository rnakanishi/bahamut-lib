#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <algorithm>

double distanceToSegment(Eigen::Array2d origin, Eigen::Array2d ending,
                         Eigen::Array2d point) {
  Eigen::Vector2d target, segment;
  target = (point - origin).matrix();
  segment = (ending - origin).matrix();
  double t = target.dot(segment) / segment.squaredNorm();
  t = std::min(1.0, std::max(0.0, t));
  Eigen::Array2d projection = origin + t * (ending - origin);

  target.normalize();
  segment.normalize();
  if (target[0] * segment[1] - target[1] * segment[0] >= 0.0)
    return (point - projection).matrix().norm();
  return -(point - projection).matrix().norm();
}
void resetLevelset(Leviathan::DualSquares &levelset) {
  auto &phi = levelset.getCellScalarData("phi");
  std::fill(phi.begin(), phi.end(), 1e8);
}

void extractLevelsetFromMesh(Leviathan::DualSquares &levelset,
                             Ramuh::LineMesh &mesh) {
  auto &phi = levelset.getCellScalarData("phi");
  auto &gradients = levelset.getCellArrayData("cellGradients");

  // For all segments
  auto &segments = mesh.getSegmentsList();
  auto &vertices = mesh.getVerticesList();
  Eigen::Vector2d tangent;
  auto h = levelset.getH();
  std::vector<bool> visited(phi.size(), false);

#pragma omp parallel for
  for (int cellId = 0; cellId < levelset.cellCount(); cellId++) {
    // Check for each cell, the closest segment
    // TOFIX: only works for convex polygons
    bool outside = false;
    for (int is = 0; is < segments.size(); is++) {
      if (!mesh.isSegmentActive(is))
        continue;
      auto segment = segments[is];

      // Get the origin point positions and find which cell it belongs
      Eigen::Array2d origin = mesh.getVertexPosition(segment[0]);
      Eigen::Array2d ending = mesh.getVertexPosition(segment[1]);
      double angle = (ending[1] - origin[1]) / (ending[0] - origin[0]);
      // auto cellId = levelset.findCellIdByCoordinate(origin);

      Eigen::Array2d cellCenter = levelset.getCellPosition(cellId);
      Eigen::Vector2d direction, target;
      direction = ending - origin;
      target = cellCenter - origin;
      double distance;

      if (cellCenter.isApprox(Eigen::Array2d(-3.7, 2.3))) {
        origin.transpose();
        // std::cerr << origin.transpose() << " " << ending.transpose()
        // << std::endl;
      }

      double cross = (target[0] * direction[1] - target[1] * direction[0]);
      int dSignal;
      if (cross == 0) {
        distance = 0.0;
        visited[cellId] = true;
        phi[cellId] = 0.0;
        gradients[cellId] = mesh.getSegmentNormal(is);
      } else {
        distance = distanceToSegment(origin, ending, cellCenter);
        if (distance > 0)
          outside = true;
        if (!visited[cellId]) {
          visited[cellId] = true;
          phi[cellId] = distance;
          gradients[cellId] = mesh.getSegmentNormal(is);
        } else {
          if (!outside) {
            phi[cellId] = std::max(phi[cellId], distance);
          } else {
            phi[cellId] = std::abs(phi[cellId]);
            if (std::abs(phi[cellId]) > std::abs(distance)) {
              phi[cellId] = std::abs(distance);
              gradients[cellId] = mesh.getSegmentNormal(is);
            }
          }
        }
      }
    }
  }
}