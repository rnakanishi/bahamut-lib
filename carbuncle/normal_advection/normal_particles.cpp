#include "normal_particles.hpp"
#include <geometry/bounding_box.h>
#include <iostream>

namespace Carbuncle {

NormalParticles3::NormalParticles3() {}

int NormalParticles3::seedParticlesOverSurface(
    Leviathan::LevelSetFluid3 levelset) {
  int particlesPerCell = 1;

  // find surface cells
  auto surfaceCells = levelset.trackSurface();
  preAllocateParticles(surfaceCells.size() * particlesPerCell);
  std::cerr << surfaceCells.size() * particlesPerCell * 2
            << " particles will be seeded\n";

  std::vector<int> normalOrigin;
  // Seed particles over the cells
  for (auto cellId : surfaceCells) {
    Ramuh::BoundingBox3 bbox = levelset.getCellBoundingBox(cellId);
    auto particles = seedParticles(bbox, particlesPerCell);

    // Attract them to the 0 levelset
    for (auto particleID : particles) {
      normalOrigin.emplace_back(particleID);
      double phi = levelset.interpolateCellScalarData(
          levelset.getCellArrayLabelId("phi"), getParticlePosition(particleID));

      double lambda = 1.0;
      int maxIt = 16;
      double threshold = 1e-6;
      Eigen::Array3d p = getParticlePosition(particleID);
      for (size_t it = 0; it < maxIt; it++) {
        auto direction = levelset.interpolateCellArrayData(
            levelset.getCellArrayLabelId("cellGradient"), p);
        phi = levelset.interpolateCellScalarData(
            levelset.getCellScalarLabelId("phi"), p);
        direction = direction.matrix().normalized().array();

        // xnew = xp + lambda(goal - phi(xp)) * N(p)
        // N(p): - levelset gradient
        Eigen::Array3d newp = p + lambda * (0.0 - phi) * direction;

        if (!levelset.getDomain().contains(newp)) {
          lambda /= 2;
          continue;
        }

        // Verify if near goal
        phi = levelset.interpolateCellScalarData(
            levelset.getCellScalarLabelId("phi"), newp);

        p = newp;
        if (std::abs(phi) < threshold)
          break;
      }
      setParticlePosition(particleID, p);
    }
  }

  for (int pid : normalOrigin) {
    auto direction = levelset.interpolateCellArrayData(
        levelset.getCellArrayLabelId("cellGradient"), getParticlePosition(pid));
    direction = direction.matrix().normalized().array();
    Eigen::Array3d normalTarget;
    normalTarget = getParticlePosition(pid) + 0.01 * direction;
    int newPid = insertParticle(normalTarget);
    normalPair[pid] = newPid;
  }

  return surfaceCells.size() * particlesPerCell * 2;
}

void NormalParticles3::print() {
  for (auto normal : normalPair) {
    auto origin = getParticlePosition(normal.first);
    auto direction = getParticlePosition(normal.second) - origin;
    std::cerr << origin.transpose() << " " << direction.transpose()
              << std::endl;
  }
}

} // namespace Carbuncle
