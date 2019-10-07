#include <fluids/particle_levelset3.h>
#include <utils/timer.hpp>
#include <blas/interpolator.h>
#include <geometry/vector3.h>
#include <cmath>
#include <set>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace Leviathan {

ParticleLevelSet3::ParticleLevelSet3()
    : ParticleLevelSet3(Eigen::Array3i(32), Ramuh::BoundingBox3::unitBox()) {}

ParticleLevelSet3::ParticleLevelSet3(Eigen::Array3i gridSize,
                                     Ramuh::BoundingBox3 domain)
    : Ramuh::ParticleSystem3(domain), Leviathan::LevelSetFluid3(gridSize,
                                                                domain) {
  _particleRadiusId = newParticleScalarLabel("radius");
  _particleVelocityId = newParticleArrayLabel("velocity");
  _particleSignalId = newParticleScalarLabel("particleSignal");
  _particleLevelSetId = newParticleScalarLabel("particleLevelSet");

  _maxParticles = 128;
}

void ParticleLevelSet3::advectEuler() {
  auto &velocity = getParticleArrayData(_particleVelocityId);
  auto &position = getParticleArrayData(_particlePositionsId);

#pragma omp parallel for
  for (size_t i = 0; i < _totalIds; i++) {
    if (_active[i]) {
      position[i] = position[i] + _dt * velocity[i];
    }
  }
}

void ParticleLevelSet3::advectParticles() {
  auto &position = getParticleArrayData(_particlePositionsId);
  std::vector<Eigen::Array3d> lastPosition(position.begin(), position.end());

  interpolateVelocityToParticles();
  advectEuler();
  interpolateVelocityToParticles();
  advectEuler();
#pragma omp parallel for
  for (size_t i = 0; i < _totalIds; i++) {
    if (_active[i])
      position[i] = 0.75 * lastPosition[i] + 0.25 * position[i];
  }
  interpolateVelocityToParticles();
  advectEuler();
#pragma omp parallel for
  for (size_t i = 0; i < _totalIds; i++) {
    if (_active[i])
      position[i] = lastPosition[i] / 3 + 2 * position[i] / 3;
  }
}

void ParticleLevelSet3::interpolateVelocityToParticles() {
  auto &velocity = ParticleSystem3::getParticleArrayData(_particleVelocityId);
  auto &position = getParticleArrayData(_particlePositionsId);
  auto h = getH();

#pragma omp parallel for
  for (size_t pid = 0; pid < _totalIds; pid++) {
    if (!_active[pid])
      continue;

    velocity[pid][0] =
        interpolateFaceScalarData(0, _faceVelocityId, position[pid]);
    velocity[pid][1] =
        interpolateFaceScalarData(1, _faceVelocityId, position[pid]);
    velocity[pid][2] =
        interpolateFaceScalarData(2, _faceVelocityId, position[pid]);
  }
}

void ParticleLevelSet3::seedParticlesNearSurface() {
  auto toSeed = findSurfaceCells(3);
  _seedCells(toSeed);
  std::cerr << "Seeded " << getParticleCount() << " particles\n";
}

void ParticleLevelSet3::_seedCells(std::set<int> &toSeed) {
  std::vector<int> cells(toSeed.begin(), toSeed.end());
  _seedCells(cells);
}

void ParticleLevelSet3::_seedCells(std::vector<int> &toSeed) {
  std::vector<int> particlesNumber(toSeed.size(), _maxParticles);
  _seedCells(toSeed, particlesNumber);
}

void ParticleLevelSet3::_seedCells(std::vector<int> &toSeed,
                                   std::vector<int> &particleNumber) {
  // For marked cells, fill them with particles
  auto &particleSignal = getParticleScalarData(_particleSignalId);
  auto &radiuses = getParticleScalarData(_particleRadiusId);
  double band[2];
  double radiusLimits[2];
  radiusLimits[0] = 0.1 * getH().maxCoeff();
  radiusLimits[1] = 0.5 * getH().maxCoeff();
  band[0] = radiusLimits[0];
  band[1] = 3.0 * getH().maxCoeff();

  for (size_t icell = 0; icell < toSeed.size(); icell++) {
    auto cell = toSeed[icell];
    auto nParticles = particleNumber[icell];

    auto box = getCellBoundingBox(cell);
    std::vector<int> seeded;

    { seeded = seedParticles(box, nParticles); }

    for (int pid : seeded) {
      Eigen::Array3d pos = getParticlePosition(pid);
      radiuses[pid] = interpolateCellScalarData(_phiId, pos);
      if ((std::abs(radiuses[pid]) < band[0]) ||
          (std::abs(radiuses[pid]) > band[1])) {
        // If seeded particle is too close from interface, remove it
        // No particle should be seeded at interface
        removeParticle(pid);
      } else {
        // Assign particle signal according to their seeded position
        if (radiuses[pid] <= 0) {
          particleSignal[pid] = -1;
        } else {
          particleSignal[pid] = 1;
        }

        radiuses[pid] =
            std::min(radiusLimits[1],
                     std::max(radiusLimits[0], std::abs(radiuses[pid])));
      }
    }
  }
} // namespace Leviathan

bool ParticleLevelSet3::_hasEscaped(int pid) {
  auto &positions = getParticleArrayData(_particlePositionsId);
  auto &signals = getParticleScalarData(_particleSignalId);
  auto &radiuses = getParticleScalarData(_particleRadiusId);
  auto h = getH();

  // FInd scaped particles, ie, particles that are on wrong side of levelset
  if (!isActive(pid))
    return false;
  Eigen::Array3d p = positions[pid];
  double levelSet = interpolateCellScalarData(_phiId, p);
  int lsignal = (levelSet <= 0) ? -1 : 1;
  int psignal = signals[pid];

  if (lsignal != psignal) {
    // Check if radius tolerance applies
    if (std::abs(levelSet) > radiuses[pid]) {
      // If their distance to interface is bigger then their radius, consider as
      // escaped
      return true;
    }
  }
  return false;
}

void ParticleLevelSet3::attractParticles() {
  auto &psignal = getParticleScalarData(_particleSignalId);
  auto &position = getParticleArrayData(_particlePositionsId);
  auto &radius = getParticleScalarData(_particleRadiusId);
  auto &particlesLevelSet = getParticleScalarData(_particleLevelSetId);

  double band[2];
  double radiusLimits[2];
  radiusLimits[0] = 0.1 * getH().maxCoeff();
  radiusLimits[1] = 0.5 * getH().maxCoeff();
  band[0] = radiusLimits[0];
  band[1] = 3.0 * getH().maxCoeff();

// Set a goal levelset for each particle
#pragma omp parallel for
  for (size_t pid = 0; pid < _totalIds; pid++) {
    if (!isActive(pid))
      continue;
    double goal;
    goal = std::abs(std::rand() % 1000000);
    goal = band[0] + goal * (band[1] - band[0]) / 1000000;
    goal = goal * psignal[pid];

    // Proceed to attraction phase
    Eigen::Array3d &p = position[pid];
    double lambda = 1.;
    double particleLevelSet = interpolateCellScalarData(_phiId, p);
    size_t maxIt = 15, it;

    for (it = 0; it < maxIt && it < maxIt; it++) {
      // Interpolate levelset
      Eigen::Array3d gradient = interpolateCellArrayData(_cellGradientId, p);

      // xnew = xp + lambda(goal - phi(xp)) * N(p)
      // N(p): - levelset gradient
      Eigen::Array3d newp = p + lambda * (goal - particleLevelSet) * gradient;

      if (!LevelSetFluid3::_domain.contains(newp)) {
        lambda /= 2;
        continue;
      }

      // Verify if near goal
      particleLevelSet = interpolateCellScalarData(_phiId, newp);
      if ((band[0] <= psignal[pid] * particleLevelSet &&
           psignal[pid] * particleLevelSet <= band[1])) {
        p = newp;
        break;
      }
      p = p + (lambda / 2) * (goal - particleLevelSet) * gradient;
    }
    if (it == maxIt) {
#pragma omp critical
      { removeParticle(pid); }
    }
  }
  adjustParticleRadius();
}

void ParticleLevelSet3::adjustParticleRadius() {
  auto &particleSignal = getParticleScalarData(_particleSignalId);
  auto &particlesLevelSet = getParticleScalarData(_particleLevelSetId);
  auto &position = getParticleArrayData(_particlePositionsId);
  auto &radius = getParticleScalarData(_particleRadiusId);
  auto &signals = getParticleScalarData(_particleSignalId);

  double radiusLimits[2];
  radiusLimits[0] = 0.1 * getH().maxCoeff();
  radiusLimits[1] = 0.5 * getH().maxCoeff();
  double limit = 3.0 * getH().maxCoeff();

  std::vector<int> toRemove;
#pragma omp parallel for
  for (size_t pid = 0; pid < _totalIds; pid++) {
    if (isActive(pid)) {
      double particleLevelSet =
          interpolateCellScalarData(_phiId, position[pid]);
      if (std::abs(particleLevelSet) > limit) {
#pragma omp critical
        { toRemove.emplace_back(pid); }
      } else {
        particlesLevelSet[pid] = particleLevelSet;
        // Set radius of each active particle according to their position in
        // relation to the interface
        if (particleLevelSet < radiusLimits[0])
          radius[pid] = radiusLimits[0];
        else if (particleLevelSet > radiusLimits[1])
          radius[pid] = radiusLimits[1];
        else
          radius[pid] = std::abs(particleLevelSet);
      }
    }
  }
  removeParticle(toRemove);
}

bool ParticleLevelSet3::correctLevelSetWithParticles() {
  auto &positions = getParticleArrayData(_particlePositionsId);
  auto &signals = getParticleScalarData(_particleSignalId);
  auto &radiuses = getParticleScalarData(_particleRadiusId);
  auto &particlesLevelSet = getParticleScalarData(_particleLevelSetId);
  auto &phi = getCellScalarData(_phiId);
  auto h = getH();

  bool hasCorrection = false;

  std::map<int, std::vector<int>> escapedCells; // contains scaped part. ids
  _escapedParticles.clear();

  // Build phi+ and phi- fields
  std::vector<double> phiPositive(phi.size()), phiNegative(phi.size());
#pragma omp parallel for
  for (size_t i = 0; i < phi.size(); i++) {
    phiPositive[i] = phiNegative[i] = phi[i];
  }

// FInd escaped particles
#pragma omp parallel for
  for (size_t pid = 0; pid < _totalIds; pid++) {
    if (!_active[pid])
      continue;
    if (_hasEscaped(pid)) {
      hasCorrection = true;
      _escapedParticles.emplace_back(pid);
      Eigen::Array3i cellIj =
          (positions[pid] - ParticleSystem3::_domain.getMin())
              .cwiseQuotient(h)
              .floor()
              .cast<int>();
      int cellId = ijkToid(cellIj[0], cellIj[1], cellIj[2]);
      escapedCells[cellId].emplace_back(pid);

      // Find which cells centers enclosure the particle
      auto particlePosition = getParticlePosition(pid);
      auto cellPosition = getCellPosition(cellId);
      if (particlePosition[0] < cellPosition[0])
        cellIj[0]--;
      if (particlePosition[1] < cellPosition[1])
        cellIj[1]--;
      if (particlePosition[2] < cellPosition[2])
        cellIj[2]--;

      std::vector<int> neighborCells;
      for (size_t dx = 0; dx < 2; dx++) {
        for (size_t dy = 0; dy < 2; dy++) {
          for (size_t dz = 0; dz < 2; dz++) {
            neighborCells.emplace_back(
                ijkToid(cellIj[0] + dx, cellIj[1] + dy, cellIj[2] + dz));
          }
        }
      }

      for (auto neighborId : neighborCells) {
        // Compute local levelset from particle
        double localPhi =
            signals[pid] *
            (radiuses[pid] - (cellPosition - particlePosition).matrix().norm());

        // Evaluate value accordign to particle sign
        if (signals[pid] > 0)
          phiPositive[cellId] = std::max(phiPositive[cellId], localPhi);
        else
          phiNegative[cellId] = std::min(phiNegative[cellId], localPhi);
      }
    }
  }
  if (hasCorrection)
#pragma omp parallel for
    for (size_t cellId = 0; cellId < cellCount(); cellId++) {
      if (std::abs(phiPositive[cellId]) < std::abs(phiNegative[cellId]))
        phi[cellId] = phiPositive[cellId];
      else
        phi[cellId] = phiNegative[cellId];
    }

  return hasCorrection;
}

int ParticleLevelSet3::findCellIdByCoordinate(Eigen::Array3d position) {
  auto h = getH();
  Eigen::Array3i cellIj = (position - ParticleSystem3::_domain.getMin())
                              .cwiseQuotient(h)
                              .floor()
                              .cast<int>();
  return ijkToid(cellIj[0], cellIj[1], cellIj[2]);
}

void ParticleLevelSet3::reseedParticles() {
  auto &signals = getParticleScalarData(_particleSignalId);
  std::vector<int> nParticles, seedingCells;

  trackSurface();
  sortParticles();
  _escapedParticles.clear();

// For all surface near cells
#pragma omp parallel for
  for (size_t imap = 0; imap < _particlesInCell.size(); imap++) {
    std::map<int, std::vector<int>>::iterator it = _particlesInCell.begin();
    std::advance(it, imap);
    int icell = it->first;
    auto particles = it->second;

    int pCount = particles.size();
    if (pCount < _maxParticles - 0.1 * _maxParticles) {
#pragma omp critical
      {
        nParticles.emplace_back(_maxParticles - pCount);
        seedingCells.emplace_back(icell);
      }
    }
    if (pCount > _maxParticles + 0.25 * _maxParticles &&
        !_isSurfaceCell[icell]) {
#pragma omp critical
      {
        removeParticle(std::vector<int>(particles.begin() + _maxParticles,
                                        particles.end()));
      }
    }
  }

  _seedCells(seedingCells, nParticles);
}

void ParticleLevelSet3::sortParticles() {
  auto h = getH();

  // Loop over particles computing their cell index
  for (size_t pid = 0; pid < getTotalParticleCount(); pid++) {
    if (!isActive(pid))
      continue;
    auto position = getParticlePosition(pid);

    // Computing cell id
    Eigen::Array3i cellIjk = (position - LevelSetFluid3::_domain.getMin())
                                 .cwiseQuotient(h)
                                 .floor()
                                 .cast<int>();
    int cellId = ijkToid(cellIjk[0], cellIjk[1], cellIjk[2]);
    if (_particlesInCell.find(cellId) == _particlesInCell.end()) {
      _particlesInCell[cellId] = std::vector<int>();
    }
    _particlesInCell[cellId].emplace_back(pid);
  }

  // for each cell that contains particles, count the particles and reseed
  for (auto particleVector : _particlesInCell) {
    auto particles = particleVector.second;
  }
}

} // namespace Leviathan
