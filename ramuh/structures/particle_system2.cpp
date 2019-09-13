#include <structures/particle_system2.h>
#include <cstdlib>

namespace Ramuh {

ParticleSystem2::ParticleSystem2() { _count = 0; }

int ParticleSystem2::particleCount() { return _count; }

int ParticleSystem2::seedParticles(BoundingBox2 region) {
  int id;
  if (_idQueue.empty()) {
    id = _active.size();
    _active.emplace_back(true);
    _positions.emplace_back(Eigen::Array2d(0));
  } else {
    id = _idQueue.front();
    _idQueue.pop();
    _active[id] = true;
  }
  Eigen::Array2d position;
  position[0] = std::rand() % 100000;
  position[1] = std::rand() % 100000;
  position = region.min() + position.cwiseProduct(region.size()) / 100000;
  _positions[id] = position;
  _count++;
  return id;
}

std::vector<int> ParticleSystem2::seedParticles(BoundingBox2 region, int n) {
  std::vector<int> ids;
  for (size_t i = 0; i < n; i++) {
    ids.emplace_back(seedParticles(region));
  }
  return ids;
}

Eigen::Array2d ParticleSystem2::getPosition(int pid) {
  if (!_active[pid])
    return Eigen::Array2d(0);
  return _positions[pid];
}

void ParticleSystem2::removeParticle(int pid) {
  if (_active[pid]) {
    _idQueue.push(pid);
    _active[pid] = false;
  }
}

void ParticleSystem2::removeParticle(std::vector<int> pids) {
  for (auto pid : pids) {
    removeParticle(pid);
  }
}

int ParticleSystem2::newScalarLabel(std::string label) {
  if (_scalarMap.find(label) == _scalarMap.end()) {
    _scalarMap[label] = _scalarData.size();
    _scalarData.emplace_back(std::vector<double>(_positions.size(), 0.));
  }
  return _scalarMap[label];
}

int ParticleSystem2::newArrayLabel(std::string label) {
  if (_arrayMap.find(label) == _arrayMap.end()) {
    _arrayMap[label] = _arrayData.size();
    _arrayData.emplace_back(
        std::vector<Eigen::Array2d>(_positions.size(), Eigen::Array2d(0.)));
  }
  return _arrayMap[label];
}

std::vector<double> &ParticleSystem2::getScalarData(std::string label) {
  return getScalarData(_scalarMap[label]);
}

std::vector<double> &ParticleSystem2::getScalarData(int id) {
  return _scalarData[id];
}

std::vector<Eigen::Array2d> &ParticleSystem2::getArrayData(std::string label) {
  return getArrayData(_arrayMap[label]);
}

std::vector<Eigen::Array2d> &ParticleSystem2::getArrayData(int id) {
  return _arrayData[id];
}

} // namespace Ramuh
