#include <structures/particle_system3.h>
#include <cstdlib>
#include <iostream>

namespace Ramuh {

ParticleSystem3::ParticleSystem3() : ParticleSystem3(BoundingBox3::unitBox()) {}

ParticleSystem3::ParticleSystem3(BoundingBox3 domain) : _domain(domain) {
  _count = 0;
  _totalIds = 0;
  _positionsId = newParticleArrayLabel("positions");

  std::srand(0);
}

int ParticleSystem3::particleCount() { return _count; }

int ParticleSystem3::seedParticles(BoundingBox3 region) {
  int id;
  if (_idQueue.empty()) {
    id = _totalIds++;
    _active.emplace_back(true);
    for (auto &dataFieldId : _scalarMap)
      _scalarData[dataFieldId.second].emplace_back(0);
    for (auto &arrayFieldId : _arrayMap)
      _arrayData[arrayFieldId.second].emplace_back(Eigen::Array3d(0));
  } else {
    id = _idQueue.front();
    _idQueue.pop();
    _active[id] = true;
  }

  auto &_positions = getParticleArrayData(_positionsId);
  Eigen::Array3d position;
  position[0] = std::rand() % 100000;
  position[1] = std::rand() % 100000;
  position[2] = std::rand() % 100000;
  position = region.min() + position.cwiseProduct(region.size()) / 100000;
  _positions[id] = position;
  _count++;
  return id;
}

std::vector<int> ParticleSystem3::seedParticles(BoundingBox3 region, int n) {
  std::vector<int> ids;
  // TODO: seedParticels: Check if domain contains seeding region
  for (size_t i = 0; i < n; i++) {
    ids.emplace_back(seedParticles(region));
  }
  return ids;
}

Eigen::Array3d ParticleSystem3::getParticlePosition(int pid) {
  if (!_active[pid])
    return Eigen::Array3d(0);
  return _arrayData[_positionsId][pid];
}

bool ParticleSystem3::isActive(int pid) { return _active[pid]; }

void ParticleSystem3::removeParticle(int pid) {
  if (_active[pid]) {
    _idQueue.push(pid);
  }
  _active[pid] = false;
  _count--;
}

void ParticleSystem3::removeParticle(std::vector<int> pids) {
  for (auto pid : pids) {
    removeParticle(pid);
  }
}

int ParticleSystem3::newParticleScalarLabel(std::string label) {
  if (_scalarMap.find(label) == _scalarMap.end()) {
    _scalarMap[label] = _scalarData.size();
    _scalarData.emplace_back(std::vector<double>(_totalIds, 0.));
  }
  return _scalarMap[label];
}

int ParticleSystem3::newParticleArrayLabel(std::string label) {
  if (_arrayMap.find(label) == _arrayMap.end()) {
    _arrayMap[label] = _arrayData.size();
    _arrayData.emplace_back(
        std::vector<Eigen::Array3d>(_totalIds, Eigen::Array3d(0.)));
  }
  return _arrayMap[label];
}

std::vector<double> &ParticleSystem3::getParticleScalarData(std::string label) {
  return getParticleScalarData(_scalarMap[label]);
}

std::vector<double> &ParticleSystem3::getParticleScalarData(int id) {
  return _scalarData[id];
}

std::vector<Eigen::Array3d> &
ParticleSystem3::getParticleArrayData(std::string label) {
  return getParticleArrayData(_arrayMap[label]);
}

std::vector<Eigen::Array3d> &ParticleSystem3::getParticleArrayData(int id) {
  return _arrayData[id];
}

std::vector<int> ParticleSystem3::searchParticles(Ramuh::BoundingBox3 region) {
  std::vector<int> particleInCell;

  for (size_t i = 0; i < _totalIds; i++) {
    if (!isActive(i))
      continue;
    Eigen::Array3d particlePosition = getParticlePosition(i);
    if (region.contains(particlePosition))
      particleInCell.emplace_back(i);
  }
  return particleInCell;
}

} // namespace Ramuh
