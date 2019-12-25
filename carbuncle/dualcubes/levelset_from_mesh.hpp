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

    // Do rasterization for each cell the segment crosses
    bool ended = false;
    Eigen::Array2d newOrigin = origin, newEnding;
    Ramuh::BoundingBox2 cellBbox = levelset.getCellBoundingBox(cellId);
    while (!ended) {
      Eigen::Vector2i cellInc(0, 0);
      newEnding = ending;
      // find which square edge intersects the surface
      if (!cellBbox.contains(ending)) {
        auto cellMin = cellBbox.getMin();
        auto cellMax = cellBbox.getMax();
        // for each edge, check intersection. If true, project to that edge:
        // bottom/top
        if (newEnding[1] < cellMin[1]) {
          newEnding[1] = cellMin[1];
          newEnding[0] = (newEnding[1] - newOrigin[1]) / angle + origin[0];
          if (cellBbox.contains(newEnding))
            cellij[1]--;
        } else if (newEnding[1] > cellMax[1]) {
          newEnding[1] = cellMax[1];
          newEnding[0] = (newEnding[1] - newOrigin[1]) / angle + origin[0];
          if (cellBbox.contains(newEnding))
            cellij[1]++;
        }
        // left/right
        if (newEnding[0] < cellMin[0]) {
          newEnding[0] = cellMin[0];
          newEnding[1] = angle * (newEnding[0] - newOrigin[0]) + newOrigin[1];
          if (cellBbox.contains(newEnding))
            cellij[0]--;
        } else if (newEnding[0] > cellMax[0]) {
          newEnding[0] = cellMax[0];
          newEnding[1] = angle * (newEnding[0] - newOrigin[0]) + newOrigin[1];
          if (cellBbox.contains(newEnding))
            cellij[0]++;
        }
      } else {
        ended = true;
      }

      // Compute cell center distance to the segment
      auto neighborCells =
          levelset.findCellNeighbors(levelset.ijToid(cellij[0], cellij[1]), 2);
      for (auto cell : neighborCells) {
        Eigen::Array2d cellCenter = levelset.getCellPosition(cell);
        Eigen::Vector2d direction, target;
        direction = ending - origin;
        target = cellCenter - origin;
        double distance;
        distance = distanceToSegment(origin, ending, cellCenter);
        if (!visited[cell]) {
          visited[cell] = true;
          phi[cell] = distance;
        } else {
          if (std::abs(phi[cell]) > std::abs(distance))
            phi[cell] = distance;
          // std::min(std::abs(phi[cell]), std::abs(distance)) * dSignal;
        }
      }
      // Special case has to be treated when segment has reach an end: check
      // for singularity
      if (ended) {
      }

      // Prepare for next cell or segment
      cellId = levelset.ijToid(cellij[0], cellij[1]);
      cellBbox = levelset.getCellBoundingBox(cellId);
      newOrigin = newEnding;
    }
  }
}