#include <structures/particle_system2.h>
#include <cstdlib>
#include <iostream>

namespace Ramuh {

ParticleSystem2::ParticleSystem2() : ParticleSystem2(BoundingBox2::unitBox()) {}

ParticleSystem2::ParticleSystem2(BoundingBox2 domain)
    : ParticleSystem2(_domain, Eigen::Array2i(22, 32)) {}

ParticleSystem2::ParticleSystem2(BoundingBox2 domain, Eigen::Array2i gridSize)
    : _domain(domain), _gridSize(gridSize) {
  _count = 0;
  _totalIds = 0;
  _particlePositionsId = newParticleArrayLabel("particlePosition");

  std::srand(0);
}

void ParticleSystem2::preAllocateParticles(int nparticles) {
  if (_totalIds < _count + nparticles) {
    _active.resize(_count + nparticles, false);
    for (auto &dataFieldId : _scalarMap)
      _scalarData[dataFieldId.second].resize(_count + nparticles, 0);
    for (auto &arrayFieldId : _arrayMap)
      _arrayData[arrayFieldId.second].resize(_count + nparticles,
                                             Eigen::Array2d(0));
  }
}

void ParticleSystem2::clearParticles() {
  for (size_t pid = 0; pid < _active.size(); pid++) {
    if (_active[pid]) {
      _active[pid] = false;
      _idQueue.push(pid);
      _count--;
    }
  }
}

int ParticleSystem2::particleCount() { return _count; }

int ParticleSystem2::insertParticle(Eigen::Array2d position) {
  int id;
  if (_idQueue.empty()) {
    id = _totalIds++;
    _active.emplace_back(true);
    for (auto &dataFieldId : _scalarMap)
      _scalarData[dataFieldId.second].emplace_back(0);
    for (auto &arrayFieldId : _arrayMap)
      _arrayData[arrayFieldId.second].emplace_back(Eigen::Array2d(0));
  } else {
    id = _idQueue.front();
    _idQueue.pop();
    _active[id] = true;
  }
  auto &_positions = getParticleArrayData(_particlePositionsId);
  _positions[id] = position;
  _count++;
  return id;
}

int ParticleSystem2::seedParticles(BoundingBox2 region) {
  int id;
  if (_idQueue.empty()) {
    id = _totalIds++;
    _active.emplace_back(true);
    for (auto &dataFieldId : _scalarMap)
      _scalarData[dataFieldId.second].emplace_back(0);
    for (auto &arrayFieldId : _arrayMap)
      _arrayData[arrayFieldId.second].emplace_back(Eigen::Array2d(0));
  } else {
    id = _idQueue.front();
    _idQueue.pop();
    _active[id] = true;
  }

  auto &_positions = getParticleArrayData(_particlePositionsId);
  Eigen::Array2d position;
  position[0] = std::rand() % 100000;
  position[1] = std::rand() % 100000;
  position = region.getMin() + position.cwiseProduct(region.getSize()) / 100000;
  _positions[id] = position;
  _count++;
  return id;
}

std::vector<int> ParticleSystem2::seedParticles(BoundingBox2 region, int n) {
  std::vector<int> ids;
  // TODO: seedParticels: Check if domain contains seeding region
  for (size_t i = 0; i < n; i++) {
    ids.emplace_back(seedParticles(region));
  }
  return ids;
}

void ParticleSystem2::setParticlePosition(int pid, Eigen::Array2d position) {
  if (_active[pid])
    _arrayData[_particlePositionsId][pid] = position;
}

Eigen::Array2d ParticleSystem2::getParticlePosition(int pid) {
  if (!_active[pid])
    return Eigen::Array2d(0);
  return _arrayData[_particlePositionsId][pid];
}

bool ParticleSystem2::isActive(int pid) { return _active[pid]; }

void ParticleSystem2::removeParticle(int pid) {
  if (_active[pid]) {
    _idQueue.push(pid);
  }
  _active[pid] = false;
  _count--;
}

void ParticleSystem2::removeParticle(std::vector<int> pids) {
  for (auto pid : pids) {
    removeParticle(pid);
  }
}

int ParticleSystem2::newParticleScalarLabel(std::string label) {
  if (_scalarMap.find(label) == _scalarMap.end()) {
    _scalarMap[label] = _scalarData.size();
    _scalarData.emplace_back(std::vector<double>(_totalIds, 0.));
  }
  return _scalarMap[label];
}

int ParticleSystem2::newParticleArrayLabel(std::string label) {
  if (_arrayMap.find(label) == _arrayMap.end()) {
    _arrayMap[label] = _arrayData.size();
    _arrayData.emplace_back(
        std::vector<Eigen::Array2d>(_totalIds, Eigen::Array2d(0.)));
  }
  return _arrayMap[label];
}

std::vector<double> &ParticleSystem2::getParticleScalarData(std::string label) {
  return getParticleScalarData(_scalarMap[label]);
}

std::vector<double> &ParticleSystem2::getParticleScalarData(int id) {
  return _scalarData[id];
}

std::vector<Eigen::Array2d> &
ParticleSystem2::getParticleArrayData(std::string label) {
  return getParticleArrayData(_arrayMap[label]);
}

std::vector<Eigen::Array2d> &ParticleSystem2::getParticleArrayData(int id) {
  return _arrayData[id];
}

std::vector<int> ParticleSystem2::searchParticles(Ramuh::BoundingBox2 region) {
  std::vector<int> particleInCell;

  for (size_t i = 0; i < _totalIds; i++) {
    if (!isActive(i))
      continue;
    Eigen::Array2d particlePosition = getParticlePosition(i);
    if (region.contains(particlePosition))
      particleInCell.emplace_back(i);
  }
  return particleInCell;
}

} // namespace Ramuh
