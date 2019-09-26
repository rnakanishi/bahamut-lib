#include <structures/particle_system2.h>
#include <cstdlib>
#include <iostream>

namespace Ramuh {

ParticleSystem2::ParticleSystem2() : ParticleSystem2(BoundingBox2::unitBox()) {}

ParticleSystem2::ParticleSystem2(BoundingBox2 domain) : _domain(domain) {
  _count = 0;
  _totalIds = 0;
  _positionsId = newParticleArrayLabel("positions");

  std::srand(0);
}

int ParticleSystem2::particleCount() { return _count; }

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

  auto &_positions = getParticleArrayData(_positionsId);
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
  // TODO: seedParticels: Check if domain contains seeding region
  for (size_t i = 0; i < n; i++) {
    ids.emplace_back(seedParticles(region));
  }
  return ids;
}

Eigen::Array2d ParticleSystem2::getParticlePosition(int pid) {
  if (!_active[pid])
    return Eigen::Array2d(0);
  return _arrayData[_positionsId][pid];
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
