#include <fluids/particle_levelset2.h>
#include <blas/interpolator.h>
#include <geometry/vector2.h>
#include <cmath>
#include <set>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace Leviathan {

ParticleLevelSet2::ParticleLevelSet2()
    : ParticleLevelSet2(Eigen::Array2i(32), Ramuh::BoundingBox2::unitBox()) {}

ParticleLevelSet2::ParticleLevelSet2(Eigen::Array2i gridSize,
                                     Ramuh::BoundingBox2 domain)
    : Ramuh::ParticleSystem2(domain), Leviathan::LevelSetFluid2(gridSize,
                                                                domain) {
  _particleRadiusId = newParticleScalarLabel("radius");
  _particleVelocityId = newParticleArrayLabel("velocity");
  _particleSignalId = newParticleScalarLabel("particleSignal");
  _particleLevelSetId = newParticleScalarLabel("particleLevelSet");
}

void ParticleLevelSet2::advectParticles() {
  auto &velocity = getParticleArrayData(_particleVelocityId);
  auto &position = getParticleArrayData(_positionsId);

  for (size_t i = 0; i < _totalIds; i++) {
    if (_active[i]) {
      position[i] = position[i] + _dt * velocity[i];
    }
    if (position[i].hasNaN() || std::isinf(position[i][0]) ||
        std::isinf(position[i][1])) {
      std::cerr << "Error";
    }
  }
}

void ParticleLevelSet2::interpolateVelocityToParticles() {
  auto &velocity = ParticleSystem2::getParticleArrayData(_particleVelocityId);
  auto &u = getFaceScalarData(0, _velocityId);
  auto &v = getFaceScalarData(1, _velocityId);
  auto &position = getParticleArrayData(_positionsId);
  auto h = getH();

  for (size_t pid = 0; pid < _totalIds; pid++) {
    if (!isActive(pid))
      continue;
    // Find which cell-id particle belongs
    Eigen::Array2i cellId = (position[pid] - ParticleSystem2::_domain.min())
                                .cwiseQuotient(h)
                                .floor()
                                .cast<int>();

    // Assemble bilinear stencil interpolation for velocities
    auto cellPos = getCellPosition(cellId[0], cellId[1]);
    int xinc, yinc;
    xinc = yinc = 1;
    if (position[pid][0] < cellPos[0])
      xinc = -1;
    if (position[pid][1] < cellPos[1])
      yinc = -1;

    std::vector<Eigen::Array2d> points;
    std::vector<double> values(4);
    double target[2];
    target[0] = position[pid][0];
    target[1] = position[pid][1];

    // u velocity
    points.emplace_back(facePosition(0, cellId[0] + 0, cellId[1] + 0));
    points.emplace_back(facePosition(0, cellId[0] + 1, cellId[1] + 0));
    points.emplace_back(facePosition(0, cellId[0] + 0, cellId[1] + yinc));
    points.emplace_back(facePosition(0, cellId[0] + 1, cellId[1] + yinc));
    values[0] = u[faceijToid(0, cellId[0] + 0, cellId[1] + 0)];
    values[1] = u[faceijToid(0, cellId[0] + 1, cellId[1] + 0)];
    values[2] = u[faceijToid(0, cellId[0] + 0, cellId[1] + yinc)];
    values[3] = u[faceijToid(0, cellId[0] + 1, cellId[1] + yinc)];
    velocity[pid][0] = Ramuh::Interpolator::bilinear(target, points, values);

    // v velocity
    points[0] = facePosition(1, cellId[0] + 0, cellId[1] + 0);
    points[1] = facePosition(1, cellId[0] + xinc, cellId[1] + 0);
    points[2] = facePosition(1, cellId[0] + 0, cellId[1] + 1);
    points[3] = facePosition(1, cellId[0] + xinc, cellId[1] + 1);
    values[0] = v[faceijToid(1, cellId[0] + 0, cellId[1] + 0)];
    values[1] = v[faceijToid(1, cellId[0] + xinc, cellId[1] + 0)];
    values[2] = v[faceijToid(1, cellId[0] + 0, cellId[1] + 1)];
    values[3] = v[faceijToid(1, cellId[0] + xinc, cellId[1] + 1)];
    velocity[pid][1] = Ramuh::Interpolator::bilinear(target, points, values);
  }
}

void ParticleLevelSet2::seedParticlesNearSurface() {
  std::vector<int> distanceToSurface(cellCount(), 1e8);
  std::vector<int> visited(cellCount(), false);
  std::set<int> toSeed;

  std::queue<int> cellQueue;
  for (auto cell : _surfaceCells) {
    distanceToSurface[cell] = 0;
    cellQueue.push(cell);
  }

  // For every cell surface tracked before, compute bfs and mark those cells
  // that are at most 3 cells away from surface
  while (!cellQueue.empty()) {
    int cell = cellQueue.front();
    cellQueue.pop();
    if (!visited[cell]) {
      visited[cell] = true;
      auto ij = idToij(cell);
      int i = ij[0], j = ij[1];
      int distance = distanceToSurface[cell];

      // All neighbors
      int neighbors[4];
      neighbors[0] = ijToid(std::max(0, i - 1), j);
      neighbors[1] = ijToid(std::min(_gridSize[0] - 1, i + 1), j);
      neighbors[2] = ijToid(i, std::max(0, j - 1));
      neighbors[3] = ijToid(i, std::min(_gridSize[1] - 1, j + 1));
      for (auto neighbor : neighbors) {
        if (!visited[neighbor]) {
          cellQueue.push(neighbor);
        }
        distanceToSurface[neighbor] =
            std::min(distanceToSurface[neighbor], distance + 1);
        if (distanceToSurface[neighbor] < 3)
          toSeed.insert(neighbor);
      }
    }
  }

  // For marked cells, fill them with particles
  auto &particleSignal = getParticleScalarData(_particleSignalId);
  int maxParticles = 128;
  for (auto cell : toSeed) {
    auto box = cellBoundingBox(cell);
    auto seeded = seedParticles(box, maxParticles);
    for (size_t i = 0; i < maxParticles / 2; i++) {
      particleSignal[seeded[i]] = -1;
    }
    for (size_t i = maxParticles / 2; i < maxParticles; i++) {
      particleSignal[seeded[i]] = +1;
    }
  }
}

void ParticleLevelSet2::attractParticles() {
  auto &psignal = getParticleScalarData(_particleSignalId);
  auto &position = getParticleArrayData(_positionsId);
  auto &radius = getParticleScalarData(_particleRadiusId);
  auto &particlesLevelSet = getParticleScalarData(_particleLevelSetId);

  // std::vector<double> goal(_totalIds, 0);
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
    Eigen::Array2d &p = position[pid];
    double lambda = 1.;
    double particleLevelSet = interpolateCellScalarData(_phiId, p);
    size_t maxIt = 15, it;

    for (it = 0; it < maxIt && it < maxIt; it++) {
      // Interpolate levelset
      Eigen::Array2d gradient = interpolateCellArrayData(_gradientId, p);

      // xnew = xp + lambda(goal - phi(xp)) * N(p)
      // N(p): - levelset gradient
      Eigen::Array2d newp = p + lambda * (goal - particleLevelSet) * gradient;

      if (!LevelSetFluid2::_domain.contains(newp)) {
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
} // namespace Leviathan

void ParticleLevelSet2::adjustParticleRadius() {
  auto &particlesLevelSet = getParticleScalarData(_particleLevelSetId);
  auto &position = getParticleArrayData(_positionsId);
  auto &radius = getParticleScalarData(_particleRadiusId);
  double radiusLimits[2];
  radiusLimits[0] = 0.1 * getH().maxCoeff();
  radiusLimits[1] = 0.5 * getH().maxCoeff();
  for (size_t pid = 0; pid < _totalIds; pid++) {
    if (isActive(pid)) {
      double particleLevelSet =
          interpolateCellScalarData(_phiId, position[pid]);
      particlesLevelSet[pid] = particleLevelSet;
      // Set radius of each active particle
      if (particleLevelSet < radiusLimits[0])
        radius[pid] = radiusLimits[0];
      else if (particleLevelSet > radiusLimits[1])
        radius[pid] = radiusLimits[1];
      else
        radius[pid] = std::abs(particleLevelSet);
    }
  }
}

void ParticleLevelSet2::correctLevelSetWithParticles() {
  auto &positions = getParticleArrayData(_positionsId);
  auto &signals = getParticleScalarData(_particleSignalId);
  auto &radiuses = getParticleScalarData(_particleRadiusId);
  auto &particleLevelSet = getParticleScalarData(_particleLevelSetId);
  auto &phi = getCellScalarData(_phiId);
  auto h = getH();

  std::vector<int> scapedParticles;
  std::map<int, std::vector<int>> scappedCells; // contains scaped part. ids
  // FInd scaped particles
  for (size_t pid = 0; pid < _totalIds; pid++) {
    if (!isActive(pid))
      continue;
    Eigen::Array2d p = positions[pid];
    double levelSet = interpolateCellScalarData(_phiId, p);
    int lsignal = (levelSet <= 0) ? -1 : 1;
    int psignal = signals[pid];

    if (lsignal != psignal) {
      // Check if radius tolerance applies
      if (std::abs(levelSet) > radiuses[pid]) {
        scapedParticles.emplace_back(pid);
        Eigen::Array2i cellIj =
            (positions[pid] - ParticleSystem2::_domain.min())
                .cwiseQuotient(h)
                .floor()
                .cast<int>();
        int cellId = ijToid(cellIj[0], cellIj[1]);
        scappedCells[cellId].emplace_back(pid);
      }
    }
  }
  {
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/particles/2d/es" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);
    auto &vel = getParticleArrayData(_particleVelocityId);
    auto &signal = getParticleScalarData(_particleSignalId);

    // for (size_t i = 0; i < particleCount(); i++) {
    for (auto i : scapedParticles) {
      if (isActive(i)) {
        auto pos = particlePosition(i);
        file << pos[0] << " " << pos[1] << " ";
        file << vel[i][0] << " " << vel[i][1] << " ";
        file << signal[i] << "\n";
      }
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }
  // Measure error and fix levelset for all scaped particles
  double phiplus, phiminus;
  for (auto cell : scappedCells) {
    int cellId = cell.first;
    auto particles = cell.second;
    Eigen::Array2d cellCenter = getCellPosition(cellId);

    phiplus = phiminus = phi[cellId];
    // For all particles in the scapped cell
    double distance[2] = {3 * h[0], 3 * h[1]}; // plus, minus
    double pradius[2] = {h[0], 0};             // plus, minus
    double finalphi[2];                        // plus, minus
    for (auto particle : particles) {
      Eigen::Array2d pPosition = particlePosition(particle);
      double pLevelset = interpolateCellScalarData(_phiId, pPosition);
      // Compute particle with least distance to center
      if (signals[particle] > 0) {
        distance[0] =
            std::min(distance[0], (pPosition - cellCenter).matrix().norm());
        pradius[0] = radiuses[particle];
      } else {
        distance[1] =
            std::min(distance[1], (pPosition - cellCenter).matrix().norm());
        pradius[1] = radiuses[particle];
      }
    }
    // Candidate levelset
    finalphi[0] = std::max(phi[cellId], (distance[0] - pradius[0]));
    finalphi[1] = std::min(phi[cellId], (-distance[1] - pradius[1]));

    // Fix phi from particles phi
    if (std::abs(finalphi[1]) < std::abs(finalphi[0]))
      phi[cellId] = finalphi[1];
    else
      phi[cellId] = finalphi[0];
  }
  adjustParticleRadius();

} // namespace Leviathan

} // namespace Leviathan
