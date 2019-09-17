#include <fluids/particle_levelset2.h>
#include <blas/interpolator.h>
#include <geometry/vector2.h>
#include <cmath>
#include <set>

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
    // Find which cell-id particle belongs
    Eigen::Array2i cellId = (position[pid] - ParticleSystem2::_domain.min())
                                .cwiseQuotient(h)
                                .floor()
                                .cast<int>();

    // Assemble bilinear stencil interpolation for velocities
    auto cellPos = cellPosition(cellId[0], cellId[1]);
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
  int maxParticles = 64;
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
} // namespace Leviathan

} // namespace Leviathan
