#include <structures/particle_system3.h>
#include <cstdlib>
#include <iostream>
#include <cmath>

namespace Ramuh {

ParticleSystem3::ParticleSystem3() : ParticleSystem3(BoundingBox3::unitBox()) {}

ParticleSystem3::ParticleSystem3(BoundingBox3 domain) : _domain(domain) {
  _count = 0;
  _totalIds = 0;
  _particlePositionsId = newParticleArrayLabel("positions");

  std::srand(0);
}

int ParticleSystem3::getParticleCount() { return _count; }

int ParticleSystem3::getTotalParticleCount() { return _totalIds; }

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

  auto &_positions = getParticleArrayData(_particlePositionsId);
  Eigen::Array3d position;
  position[0] = std::rand() % 100000;
  position[1] = std::rand() % 100000;
  position[2] = std::rand() % 100000;
  position = region.getMin() + position.cwiseProduct(region.getSize()) / 100000;
  _positions[id] = position;
  _count++;
  return id;
}

void ParticleSystem3::preAllocateParticles(int nparticles) {
  if (_totalIds < _count + nparticles) {
    _active.resize(_count + nparticles, false);
    for (auto &dataFieldId : _scalarMap)
      _scalarData[dataFieldId.second].resize(_count + nparticles, 0);
    for (auto &arrayFieldId : _arrayMap)
      _arrayData[arrayFieldId.second].resize(_count + nparticles,
                                             Eigen::Array3d(0));
  }
}

std::vector<int> ParticleSystem3::seedParticles(BoundingBox3 region, int n) {
  auto &_positions = getParticleArrayData(_particlePositionsId);
  std::vector<int> ids(n);

  if (_count + n >= _scalarData[0].size()) {
    std::cerr << "\033[21;33m[WARNING]: \033[0mData not preallocated!\n";
  }

  for (size_t i = 0; i < n; i++) {
    Eigen::Array3d position;
#pragma omp critical(insertion)
    {
      if (_idQueue.empty()) {
        ids[i] = _totalIds++;
      } else {
        ids[i] = _idQueue.front();
        _idQueue.pop();
      }
      _active[ids[i]] = true;

      position[0] = std::rand() % 100000;
      position[1] = std::rand() % 100000;
      position[2] = std::rand() % 100000;
      position =
          region.getMin() + position.cwiseProduct(region.getSize()) / 100000;
      _positions[ids[i]] = position;
    }
#pragma omp atomic
    _count++;
  }

  return ids;
}

Eigen::Array3d ParticleSystem3::getParticlePosition(int pid) {
  if (!_active[pid])
    return Eigen::Array3d(0);
  return _arrayData[_particlePositionsId][pid];
}

bool ParticleSystem3::isActive(int pid) { return _active[pid]; }

void ParticleSystem3::removeParticle(int pid) {
  if (_active[pid]) {
    _idQueue.push(pid);
    _count--;
  }
  _active[pid] = false;
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
