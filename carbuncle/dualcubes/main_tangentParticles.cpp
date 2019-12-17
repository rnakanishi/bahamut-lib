#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <iomanip>
#include "dual_cubes.h"
#include <surface/dual_squares.h>
#include <structures/line_mesh.hpp>
#include "initialization.hpp"
#include "tangent_particles2.hpp"

void printParticles(Carbuncle::TangentParticles2 particles) {
  std::ofstream file("matlab/particles.txt");
  auto &normalPair = particles.getPairMap();
  for (auto normal : normalPair) {
    if (particles.isActive(normal.first) && particles.isActive(normal.second)) {
      auto origin = particles.getParticlePosition(normal.first);
      auto ending = particles.getParticlePosition(normal.second);
      Eigen::Array2d tangent((ending - origin).matrix());
      Eigen::Vector2d direction(tangent[1], -tangent[0]);
      file << origin.transpose() << " " << direction.normalized().transpose()
           << std::endl;
    }
  }
  file.close();
}

double distanceToSegment(Eigen::Array2d origin, Eigen::Array2d ending,
                         Eigen::Array2d point) {
  double length = (origin - point).matrix().squaredNorm();
  Eigen::Vector2d target, segment;
  target = (point - origin).matrix();
  segment = (ending - origin).matrix();
  double t = target.dot(segment) / length;
  t = std::min(1.0, std::max(0.0, t));
  Eigen::Array2d projection = origin + t * (ending - origin);

  target.normalize();
  segment.normalize();
  if (target[0] * segment[1] - target[1] * segment[0] > 0.0)
    return (point - projection).matrix().norm();
  return -(point - projection).matrix().norm();
}

void fixLevelsetFromMesh(Leviathan::DualSquares &levelset,
                         Ramuh::LineMesh &mesh) {
  auto &phi = levelset.getCellScalarData("phi");

  // For all segments
  auto &segments = mesh.getSegmentsList();
  auto &vertices = mesh.getVerticesList();
  Eigen::Vector2d tangent;
  auto h = levelset.getH();

  for (auto segment : segments) {
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
    std::vector<bool> visited(phi.size(), false);
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
      auto neighborCells = levelset.findCellNeighbors(cellId, 2);
      for (auto cell : neighborCells) {
        Eigen::Array2d cellCenter = levelset.getCellPosition(cell);
        Eigen::Vector2d direction, target;
        direction = newEnding - newOrigin;
        target = cellCenter - newOrigin;
        int pSignal = (phi[cell] < 0) ? -1 : 1;
        double distance;
        if ((target[0] * direction[1] - target[1] * direction[0]) * pSignal < 0)
          distance = distanceToSegment(newOrigin, newEnding, cellCenter);
        else
          distance = distanceToSegment(newEnding, newOrigin, cellCenter);
        int dSignal = (distance < 0) ? -1 : 1;
        // if (std::abs(phi[cell] > std::abs(distance)))
        if (dSignal == pSignal)
          if (!visited[cell]) {
            visited[cell] = true;
            phi[cell] = distance;
          } else {
            phi[cell] =
                std::min(std::abs(phi[cell]), std::abs(distance)) * pSignal;
          }
      }
      // Special case has to be treated when segment has reach an end: check for
      // singularity
      if (ended) {
      }

      // Prepare for next cell or segment
      cellId = levelset.ijToid(cellij[0], cellij[1]);
      cellBbox = levelset.getCellBoundingBox(cellId);
      newOrigin = newEnding;
    }
  }
}

int main(int argc, char const *argv[]) {
  Leviathan::DualSquares cubes(
      Eigen::Array2i(40, 40),
      Ramuh::BoundingBox2(Eigen::Array2d(-5, -5), Eigen::Array2d(5, 5)));
  Carbuncle::TangentParticles2 particles(cubes);

  Eigen::Array2d center = Eigen::Array2d(0, 0);
  double radius = 1.45;

  initializeCube(cubes, center, radius, ParametricSurface::SQUARE);
  initializeGradientsAtIntersection(cubes, center, radius,
                                    ParametricSurface::SQUARE);
  cubes.setFolder("results/dualSquares/");
  Ramuh::LineMesh squareMesh;
  createLineMesh(squareMesh, center, radius);

  defineCellsVelocity(cubes);

  particles.seedParticlesOverSurface(cubes, squareMesh);
  defineParticlesVelocity(particles);
  printParticles(particles);

  // cubes.computeCellsGradient();
  // cubes.computeIntersectionAndNormals();

  cubes.extractSurface();
  squareMesh = particles.extractSurface(cubes);

  for (int i = 1; i <= 50; i++) {
    // particles.clearParticles();
    // particles.seedParticlesOverSurface(cubes, squareMesh);

    // fixLevelsetFromMesh(cubes, squareMesh);
    // cubes.redistance();
    // cubes.computeCellsGradient();
    cubes.advectWeno();
    particles.advectParticles();

    // After levelset and particles advection, cell gradient should be corrected
    // using particle information
    // Correct computed gradients with particle normals
    particles.fixLevelsetGradients(cubes);
    // cubes.computeCellsGradient();
    cubes.computeIntersectionAndNormals();
    // particles.estimateCellNormals(cubes);
    cubes.extractSurface();
    squareMesh = particles.extractSurface(cubes);
    printParticles(particles);
  }
  // cubes.computeCellsGradient();
  // cubes.computeIntersectionAndNormals();
  // particles.fixLevelsetGradients(cubes);
  // cubes.extractSurface();

  particles.extractSurface(cubes);
  return 0;
}
