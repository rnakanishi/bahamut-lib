#include "normal_particles2.hpp"
#include <geometry/bounding_box.h>
#include <geometry/dual_marching.h>
#include <blas/interpolator.h>
#include <iostream>

namespace Carbuncle {

NormalParticles2::NormalParticles2() {}

NormalParticles2::NormalParticles2(Leviathan::LevelSetFluid2 levelset)
    : ParticleLevelSet2(levelset.getGridSize(), levelset.getDomain()) {}

std::map<int, int> &NormalParticles2::getPairMap() { return normalPair; }

int NormalParticles2::seedParticlesOverSurface(
    Leviathan::LevelSetFluid2 levelset) {
  // TODO: Seed particles only near the cell face
  // TODO: Also ensure uniform sampling of the particles along surface
  // Check which cell it intersect with and seed particles there
  int particlesPerCell = 1;
  auto gradient = levelset.getCellArrayData("cellGradient");
  // find surface cells
  auto surfaceCells = levelset.trackSurface();
  // preAllocateParticles(surfaceCells.size() * particlesPerCell);
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
      double phi = levelset.interpolateCellScalarData(
          "phi", getParticlePosition(particleID));

      double lambda = 1.0;
      int maxIt = 16;
      double threshold = 1e-6;
      Eigen::Array2d p = getParticlePosition(particleID);
      for (size_t it = 0; it < maxIt; it++) {
        // auto direction = levelset.interpolateCellArrayData("cellGradient",
        // p);
        auto direction = gradient[cellId];
        phi = levelset.interpolateCellScalarData("phi", p);
        direction = direction.matrix().normalized().array();

        // xnew = xp + lambda(goal - phi(xp)) * N(p)
        // N(p): - levelset gradient
        Eigen::Array2d newp = p + lambda * (0.0 - phi) * direction;

        if (!levelset.getDomain().contains(newp)) {
          lambda /= 2;
          continue;
        }

        // Verify if near goal
        phi = levelset.interpolateCellScalarData("phi", newp);

        p = newp;
        if (std::abs(phi) < threshold)
          break;
      }
      setParticlePosition(particleID, p);

      // For each added particle, add a pair for normal target
      // copy the cell gradient to the particle
      auto direction = gradient[cellId];
      direction = direction.matrix().normalized().array();
      Eigen::Array2d normalTarget;
      normalTarget = getParticlePosition(particleID) + 0.01 * direction;
      int newPid = insertParticle(normalTarget);
      normalPair[particleID] = newPid;
    }
  }

  return surfaceCells.size() * particlesPerCell * 2;
}

void NormalParticles2::estimateCellNormals(
    Leviathan::LevelSetFluid2 &levelset) {
  auto h = levelset.getH();
  auto &gradients = levelset.getCellArrayData("cellGradient");
  std::map<int, std::vector<int>> cellMap; // store particles id in the cell

  // Iterate through all particles and find their corresponding cell
  for (auto &particle : normalPair) {
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
      Eigen::Array2d ending = getParticlePosition(normalPair[pid]);

      normalPositions.emplace_back(origin);
      normalVectors.emplace_back((ending - origin).matrix());
    }

    // Interpolate the final vector from the particles vector
    Eigen::Vector2d finalNormal;
    auto cellPosition = levelset.getCellPosition(cell.first);
    finalNormal = Ramuh::Interpolator::shepard(cellPosition, normalPositions,
                                               normalVectors);
    gradients[cell.first] = finalNormal.normalized().array();
  }
}

void NormalParticles2::advectParticles() {
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

void NormalParticles2::defineParticlesVelocity() {
  auto &velocity = getParticleArrayData("particleVelocity");
  for (size_t pid = 0; pid < particleCount(); pid++) {
    auto p = getParticlePosition(pid);
    velocity[pid] = Eigen::Array2d(-p[1], p[0]);
  }
}

void NormalParticles2::fixLevelsetGradients(
    Leviathan::LevelSetFluid2 &levelset) {
  auto &gradient = levelset.getCellArrayData("cellGradient");
  auto &unormalPosition = levelset.getFaceArrayData(0, "facePosition");
  auto &vnormalPosition = levelset.getFaceArrayData(1, "facePosition");
  auto &unormal = levelset.getFaceArrayData(0, "faceNormal");
  auto &vnormal = levelset.getFaceArrayData(1, "faceNormal");

  // Track all surface cells
  auto surfaceCells = levelset.trackSurface();

  // Find which cells contain particles
  // Interpolate gradient values from the closest particles

  // For those cells that may not contain any particle, copy from closest
  // particle
}

void NormalParticles2::extractSurface(Leviathan::LevelSetFluid2 &levelset) {

  // Find surface cells
  auto surfaceCells = levelset.trackSurface();
  std::vector<Eigen::Array2d> positions;
  std::vector<Eigen::Vector2d> normals;

  Ramuh::DualMarching2 surface;

  surface.setBaseFolder("results/dualSquares/particles/");

  // Check which particles belong to that cell
  sortParticles();
  surfaceCells = computeCellsWithParticles();
  int nParticles = particleCount() / 2;
  for (auto cell : surfaceCells) {
    auto particles = getParticlesInCell(cell);

    positions.clear();
    normals.clear();
    // Use particles and their normals to compute minimization point for surface
    for (auto particle : particles) {
      if (normalPair.find(particle) == normalPair.end())
        continue;
      Eigen::Vector2d normal;
      Eigen::Array2d position = getParticlePosition(particle);
      Eigen::Array2d target = getParticlePosition(normalPair[particle]);
      normal = target - position;

      positions.emplace_back(position);
      normals.emplace_back(normal);
    }
    if (positions.size() == 1)
      surface.getPoints().emplace_back(positions[0]);
    else if (!positions.empty()) {
      auto ij = levelset.idToij(cell);
      Eigen::Array2i index(ij[0], ij[1]);
      surface.evaluateSquare(index, positions, normals,
                             levelset.getCellBoundingBox(cell));
    }
  }
  surface.reconstruct();
}

} // namespace Carbuncle
