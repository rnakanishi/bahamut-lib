#ifndef __RAMUH_PARTICLE_SYSTEM_3_H__
#define __RAMUH_PARTICLE_SYSTEM_3_H__

#include <Eigen/Dense>
#include <geometry/bounding_box.h>
#include <vector>
#include <queue>
#include <map>
#include <string>

namespace Ramuh {

class ParticleSystem3 {
public:
  ParticleSystem3();
  ParticleSystem3(BoundingBox3 domain);
  ParticleSystem3(BoundingBox3 domain, Eigen::Array3i gridSize);

  int getParticleCount();
  int getTotalParticleCount();

  /**
   * @brief Seeds a predertmined number (default: 1) of particles over a given
   * region and set them to active. If there are non-active ids, they are used
   * first. The spread of all particles is randomly generated.
   *
   * @param region bounding box corresponding to target region
   * @param n [optional] number of particles to be seeded
   * @return int Particle ids
   */
  int seedParticles(BoundingBox3 region);
  std::vector<int> seedParticles(BoundingBox3 region, int n);

  /**
   * @brief Get the Positions of a single particle. If the required particle is
   * not active, then origin coordinate is returned.
   *
   * @param pid id of the particle
   * @return std::vector<Eigen::Array3d> Either a vector containing all particle
   * positions; or a single position
   */
  Eigen::Array3d getParticlePosition(int pid);

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

  std::vector<int> searchParticles(Ramuh::BoundingBox3 region);

  int newParticleScalarLabel(std::string label);

  int newParticleArrayLabel(std::string label);

  std::vector<double> &getParticleScalarData(std::string label);
  std::vector<double> &getParticleScalarData(int id);

  std::vector<Eigen::Array3d> &getParticleArrayData(std::string label);
  std::vector<Eigen::Array3d> &getParticleArrayData(int id);

protected:
  std::vector<bool> _active;
  std::queue<int> _idQueue;

  int _particlePositionsId;

  std::vector<std::vector<double>> _scalarData;
  std::vector<std::vector<Eigen::Array3d>> _arrayData;
  std::map<std::string, int> _arrayMap, _scalarMap;

  BoundingBox3 _domain;
  Eigen::Array3i _gridSize;
  int _count, _totalIds;
};

} // namespace Ramuh

#endif