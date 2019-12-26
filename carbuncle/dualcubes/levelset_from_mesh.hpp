#include <geometry/bounding_box.h>
#include <Eigen/Dense>

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
    for (int is = 0; is < segments.size(); is++) {
      if (!mesh.isSegmentActive(is))
        continue;
      auto segment = segments[is];

      // Get the origin point positions and find which cell it belongs
      Eigen::Array2d origin = mesh.getVertexPosition(segment[0]);
      Eigen::Array2d ending = mesh.getVertexPosition(segment[1]);
      double angle = (ending[1] - origin[1]) / (ending[0] - origin[0]);
      auto cellId = levelset.findCellIdByCoordinate(origin);
      auto cellij = levelset.idToij(cellId);

      Eigen::Array2d cellCenter = levelset.getCellPosition(cellId);
      Eigen::Vector2d direction, target;
      direction = ending - origin;
      target = cellCenter - origin;
      double distance;
      distance = distanceToSegment(origin, ending, cellCenter);
      if (!visited[cellId]) {
        visited[cellId] = true;
        phi[cellId] = distance;
      } else {
        if (std::abs(phi[cellId]) > std::abs(distance))
          phi[cellId] = distance;
        // std::min(std::abs(phi[cell]), std::abs(distance)) * dSignal;
      }
    }
  }
}