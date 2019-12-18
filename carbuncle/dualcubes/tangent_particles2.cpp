#include "tangent_particles2.hpp"
#include <geometry/bounding_box.h>
#include <geometry/dual_marching.h>
#include <blas/interpolator.h>
#include <iostream>
#include <cstdlib>

namespace Carbuncle {

double getFunctionValue(Eigen::Array2d p) {
  Eigen::Array2d center(0, 0);
  double radius = 1.5;

  double x, y, x2, y2;
  double distance;
  x = p[0] - center[0];
  y = p[1] - center[1];
  x2 = x * x;
  y2 = y * y;

  // CUBE
  distance = std::max(std::fabs(x), std::fabs(y)) - radius;
  if (distance > 0) {
    p = p.abs();
    distance = 0.0;
    x = std::max(0.0, p[0] - radius);
    y = std::max(0.0, p[1] - radius);
    distance = sqrt(x * x + y * y);
  }
  return distance;
}

TangentParticles2::TangentParticles2() {}

TangentParticles2::TangentParticles2(Leviathan::LevelSetFluid2 levelset)
    : ParticleLevelSet2(levelset.getGridSize(), levelset.getDomain()) {}

std::map<int, int> &TangentParticles2::getPairMap() { return tangentPair; }

int TangentParticles2::seedParticlesOverSurface(
    Leviathan::LevelSetFluid2 levelset) {
  // TODO: Seed particles only near the cell face
  // TODO: Also ensure uniform sampling of the particles along surface
  // Check which cell it intersect with and seed particles there
  int particlesPerCell = 15;
  auto gradient = levelset.getCellArrayData("cellGradient");
  // find surface cells
  auto surfaceCells = levelset.trackSurface();
  std::cerr << surfaceCells.size() * particlesPerCell * 2
            << " particles will be seeded\n";

  std::vector<int> normalOrigin;
  // Seed particles over the cells
  for (auto cellId : surfaceCells) {
    Ramuh::BoundingBox2 bbox = levelset.getCellBoundingBox(cellId);
    auto particles = seedParticles(bbox, particlesPerCell);

    // Attract them to the 0 levelset
    for (auto particleID : particles) {
      normalOrigin.emplace_back(particleID);
      double phi;

      double lambda = 1.0;
      int maxIt = 16;
      double threshold = 1e-6;
      Eigen::Array2d p = getParticlePosition(particleID);
      for (size_t it = 0; it < maxIt; it++) {
        auto direction = gradient[cellId];

        // Analytic value is being used to get phi value
        phi = getFunctionValue(p);
        direction = direction.matrix().normalized().array();

        Eigen::Array2d newp = p + lambda * (0.0 - phi) * direction;
        if (!levelset.getDomain().contains(newp)) {
          // In case the new point falls of the domain
          lambda /= 2;
          continue;
        }

        // Verify if near goal
        p = newp;
        phi = getFunctionValue(newp);
        if (std::abs(phi) < threshold)
          break;
      }
      setParticlePosition(particleID, p);

      // For each added particle, add a pair for tangent vector with simplified
      // cross product: (vx, vy, 0) x (0, 0, 1)
      auto direction = gradient[cellId];
      Eigen::Array2d tangent;
      direction = direction.matrix().normalized().array();
      tangent[0] = direction[1];
      tangent[1] = -direction[0];
      Eigen::Array2d target;
      target = getParticlePosition(particleID) + 0.01 * tangent;
      int newPid = insertParticle(target);
      tangentPair[particleID] = newPid;
    }
  }

  return surfaceCells.size() * particlesPerCell * 2;
}

int TangentParticles2::seedParticlesOverSurface(
    Leviathan::LevelSetFluid2 levelset, Ramuh::LineMesh mesh) {
  // For all segments
  auto &segments = mesh.getSegmentsList();
  auto &vertices = mesh.getVerticesList();
  Eigen::Vector2d tangent;
  auto h = levelset.getH();
  int nTanPoints = 25;

  for (auto segment : segments) {
    // Get the origin point positions and find which cell it belongs
    Eigen::Array2d origin = mesh.getVertexPosition(segment[0]);
    Eigen::Array2d ending = mesh.getVertexPosition(segment[1]);
    tangent = (ending - origin).matrix().normalized();
    double angle = (ending[1] - origin[1]) / (ending[0] - origin[0]);
    auto cellId = levelset.findCellIdByCoordinate(origin);
    auto cellij = levelset.idToij(cellId);

    // Do rasterization for each cell the segment crosses
    // Generate one tangetnPoint for each cell-segment pair
    bool ended = false;
    Eigen::Array2d newOrigin = origin, newEnding;
    Ramuh::BoundingBox2 cellBbox = levelset.getCellBoundingBox(cellId);
    while (!ended) {
      Eigen::Vector2i cellInc(0, 0);
      newEnding = ending;
      if (!cellBbox.contains(ending)) {
        // find which square edge intersects the surface
        auto cellMin = cellBbox.getMin();
        auto cellMax = cellBbox.getMax();
        // for each edge, check intersection. If true, project to that edge
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
        cellId = levelset.ijToid(cellij[0], cellij[1]);
        cellBbox = levelset.getCellBoundingBox(cellId);
      } else {
        ended = true;
      }

      // Generate a certain quantity of point along the cell segment
      Eigen::Array2d newPoint, tangentTarget;
      tangent = (newEnding - newOrigin).matrix();
      Eigen::Array2d samplOrigin = newOrigin + tangent.array() * 0.001;
      Eigen::Array2d samplEnding = newEnding - tangent.array() * 0.001;
      Eigen::Vector2d segTangent = (samplEnding - samplOrigin).matrix();
      for (int samplI = 0; samplI <= nTanPoints; samplI++) {
        newPoint = samplOrigin +
                   segTangent.array() * (1.0 - (double)samplI / nTanPoints);
        tangentTarget = newPoint + (segTangent.array() * 0.001);
        int pid = insertParticle(newPoint);
        int tid = insertParticle(tangentTarget);
        tangentPair[pid] = tid;
      }

      newOrigin = newEnding;
    }
  }
  return particleCount();
}

void TangentParticles2::estimateCellNormals(
    Leviathan::LevelSetFluid2 &levelset) {
  auto h = levelset.getH();
  auto &gradients = levelset.getCellArrayData("cellGradient");
  std::map<int, std::vector<int>> cellMap; // store particles id in each cell

  // Iterate through all particles and find their corresponding cell
  for (auto &particle : tangentPair) {
    // Only need to do it for the origin of the vector
    int pid = particle.first;
    auto position = getParticlePosition(pid);
    Eigen::Array2i cellIj = (position - levelset.getDomain().getMin())
                                .cwiseQuotient(h)
                                .floor()
                                .cast<int>();
    int cellId = levelset.ijToid(cellIj[0], cellIj[1]);

    if (cellMap.find(cellId) == cellMap.end())
      cellMap[cellId] = std::vector<int>();
    cellMap[cellId].emplace_back(pid);
  }

  // For all cells that contains particles:
  for (auto cell : cellMap) {
    // Get the cell center position
    auto cellCenter = levelset.getCellPosition(cell.first);

    // Compute all vectors for that cell
    std::vector<Eigen::Vector2d> normalVectors;
    std::vector<Eigen::Array2d> normalPositions;
    for (auto pid : cell.second) {
      Eigen::Array2d origin = getParticlePosition(pid);
      Eigen::Array2d ending = getParticlePosition(tangentPair[pid]);
      Eigen::Array2d tangent((ending - origin).matrix());

      normalPositions.emplace_back(origin);
      normalVectors.emplace_back(tangent[1], -tangent[0]);
    }

    // Interpolate the final vector from the particles vector
    Eigen::Vector2d finalNormal;
    auto cellPosition = levelset.getCellPosition(cell.first);
    finalNormal = Ramuh::Interpolator::closestPoint(
        cellPosition, normalPositions, normalVectors);
    gradients[cell.first] = finalNormal.normalized().array();
  }
}

void TangentParticles2::advectParticles() {
  auto &position = getParticleArrayData(_particlePositionsId);

  std::vector<Eigen::Array2d> lastPosition(position.size());
  for (size_t i = 0; i < position.size(); i++) {
    lastPosition[i] = position[i];
  }
  advectEuler();
  defineParticlesVelocity();
  advectEuler();
#pragma omp parallel for
  for (size_t i = 0; i < position.size(); i++) {
    position[i] = 0.75 * lastPosition[i] + 0.25 * position[i];
  }
  defineParticlesVelocity();
  advectEuler();
#pragma omp parallel for
  for (size_t i = 0; i < position.size(); i++) {
    position[i] = lastPosition[i] / 3 + 2 * position[i] / 3;
  }
}

void TangentParticles2::defineParticlesVelocity() {
  auto &velocity = getParticleArrayData("particleVelocity");
  for (size_t pid = 0; pid < particleCount(); pid++) {
    auto p = getParticlePosition(pid);
    velocity[pid] = Eigen::Array2d(-p[1], p[0]);
  }
}

void TangentParticles2::fixLevelsetGradients(
    Leviathan::LevelSetFluid2 &levelset) {
  auto &gradient = levelset.getCellArrayData("cellGradient");
  auto &unormalPosition = levelset.getFaceArrayData(0, "facePosition");
  auto &vnormalPosition = levelset.getFaceArrayData(1, "facePosition");
  auto &unormal = levelset.getFaceArrayData(0, "faceNormal");
  auto &vnormal = levelset.getFaceArrayData(1, "faceNormal");

  // Track all surface cells
  auto surfaceCells = levelset.trackSurface();
  auto h = levelset.getH();

  std::vector<Eigen::Array2d> allPositions;
  std::vector<Eigen::Vector2d> allNormals;

  for (size_t pid = 0; pid < particleCount(); pid++) {
    if (tangentPair.find(pid) != tangentPair.end()) {
      Eigen::Array2d origin = getParticlePosition(pid);
      Eigen::Array2d ending = getParticlePosition(tangentPair[pid]);
      Eigen::Array2d tangent((ending - origin).matrix());

      allPositions.emplace_back(origin);
      allNormals.emplace_back(-tangent[1], tangent[0]);
    }
  }
  if (allPositions.empty())
    return;

  // for all cells that is surface, fix its face normals
  for (auto cell : surfaceCells) {
    auto ij = levelset.idToij(cell);

    // Closest point normal is used since interpolation may smooth the vectors
    unormal[levelset.faceijToid(0, ij[0], ij[1])] =
        Ramuh::Interpolator::closestPoint(
            unormalPosition[levelset.faceijToid(0, ij[0], ij[1])], allPositions,
            allNormals)
            .array();
    vnormal[levelset.faceijToid(1, ij[0], ij[1])] =
        Ramuh::Interpolator::closestPoint(
            vnormalPosition[levelset.faceijToid(1, ij[0], ij[1])], allPositions,
            allNormals)
            .array();
  }
}

Ramuh::LineMesh
TangentParticles2::extractSurface(Leviathan::LevelSetFluid2 &levelset) {
  if (particleCount() == 0)
    return Ramuh::LineMesh();

  // Find surface cells
  auto surfaceCells = levelset.trackSurface();
  std::vector<Eigen::Array2d> positions;
  std::vector<Eigen::Vector2d> normals;
  static Ramuh::DualMarching2 surface(levelset.getResolution());

  surface.setBaseFolder("results/dualSquares/particles/");
  surface.clear();

  // Check which particles belong to that cell
  surfaceCells = computeCellsWithParticles();

  for (auto cell : surfaceCells) {
    auto particles = getParticlesInCell(cell);
    positions.clear();
    normals.clear();
    // Use particles and their normals to compute minimization point for surface
    for (auto particle : particles) {
      if (tangentPair.find(particle) == tangentPair.end())
        continue;
      Eigen::Array2d position = getParticlePosition(particle);
      Eigen::Array2d target = getParticlePosition(tangentPair[particle]);
      Eigen::Array2d tangent((target - position).matrix());
      Eigen::Vector2d normal(tangent[1], -tangent[0]);

      positions.emplace_back(position);
      normals.emplace_back(normal.normalized());
    }
    if (positions.size() < 1) {
      surface.getPoints().emplace_back(positions[0]);
    } else if (!positions.empty()) {
      auto ij = levelset.idToij(cell);
      Eigen::Array2i index(ij[0], ij[1]);
      surface.evaluateSquare(index, positions, normals);
    }
    // previousCell = cell;
  }
  return surface.reconstruct();
}

} // namespace Carbuncle
