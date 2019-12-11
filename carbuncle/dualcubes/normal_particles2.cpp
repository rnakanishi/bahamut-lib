#include "normal_particles2.hpp"
#include <geometry/bounding_box.h>
#include <blas/interpolator.h>
#include <iostream>

namespace Carbuncle {

NormalParticles2::NormalParticles2() {}

std::map<int, int> &NormalParticles2::getPairMap() { return normalPair; }

int NormalParticles2::seedParticlesOverSurface(
    Leviathan::LevelSetFluid2 levelset) {
  int particlesPerCell = 1;

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
        auto direction = levelset.interpolateCellArrayData("cellGradient", p);
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
    }
  }

  for (int pid : normalOrigin) {
    auto direction = levelset.interpolateCellArrayData(
        "cellGradient", getParticlePosition(pid));
    direction = direction.matrix().normalized().array();
    Eigen::Array2d normalTarget;
    normalTarget = getParticlePosition(pid) + 0.01 * direction;
    int newPid = insertParticle(normalTarget);
    normalPair[pid] = newPid;
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

} // namespace Carbuncle
