#ifndef __RAMUH_PARTICLE_SYSTEM_2_HPP
#define __RAMUH_PARTICLE_SYSTEM_2_HPP

#include <Eigen/Dense>
#include <geometry/bounding_box.h>
#include <vector>
#include <queue>
#include <map>
#include <string>

namespace Ramuh {

class ParticleSystem2 {
public:
  ParticleSystem2();
  ParticleSystem2(BoundingBox2 domain);
  ParticleSystem2(BoundingBox2 domain, Eigen::Array2i gridSize);

  int particleCount();

  /**
   * @brief Seeds a predertmined number (default: 1) of particles over a given
   * region and set them to active. If there are non-active ids, they are used
   * first. The spread of all particles is randomly generated.
   *
   * @param region bounding box corresponding to target region
   * @param n [optional] number of particles to be seeded
   * @return int Particle ids
   */
  int seedParticles(BoundingBox2 region);
  std::vector<int> seedParticles(BoundingBox2 region, int n);

  /**
   * @brief Get the Positions of a single particle. If the required particle is
   * not active, then origin coordinate is returned.
   *
   * @param pid id of the particle
   * @return std::vector<Eigen::Array2d> Either a vector containing all particle
   * positions; or a single position
   */
  Eigen::Array2d getParticlePosition(int pid);

  /**
   * @brief Return either a particle is active or not.
   *
   * @param pid particle id
   * @return true if particle is active
   * @return false if particle was removed from the structure
   */
  bool isActive(int pid);

  /**
   * @brief Remove particles from the particle system. Parameter can be eigher a
   * single id, or a vector with id bundle. Instead of actually removing them
   * from the data vector, only  mark them as inactive and theur ids are queued
   * to be used when new particles are inserted.
   *
   * @param pid either a single particle id, or a vector with id bundle
   */
  void removeParticle(int pid);
  void removeParticle(std::vector<int> pids);

  std::vector<int> searchParticles(Ramuh::BoundingBox2 region);

  int newParticleScalarLabel(std::string label);

  int newParticleArrayLabel(std::string label);

  std::vector<double> &getParticleScalarData(std::string label);
  std::vector<double> &getParticleScalarData(int id);

  std::vector<Eigen::Array2d> &getParticleArrayData(std::string label);
  std::vector<Eigen::Array2d> &getParticleArrayData(int id);

protected:
  std::vector<bool> _active;
  std::queue<int> _idQueue;
  int _positionsId;

  std::vector<std::vector<double>> _scalarData;
  std::vector<std::vector<Eigen::Array2d>> _arrayData;
  std::map<std::string, int> _arrayMap, _scalarMap;

  Eigen::Array2i _gridSize;

  BoundingBox2 _domain;
  int _count, _totalIds;
};

} // namespace Ramuh
#endif