#ifndef __LEVIATHAN_LEVELSET_FLUID_H__
#define __LEVIATHAN_LEVELSET_FLUID_H__

#include <structures/mac_grid3.h>
#include <fstream>
#include <sstream>

namespace Leviathan {

class LevelSetFluid3 : public Ramuh::MacGrid3 {
public:
  LevelSetFluid3();
  LevelSetFluid3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain);

  /**
   * @brief compute semi Lagrangean advection for levelset values.
   *
   */
  void advectSemiLagrangean();
  void advect(int order);

  void advectMacCormack();

  void advectWeno();

  void advectUpwind();

  void computeCellsGradient();

  void redistance();

protected:
  double __interpolateVelocityU(Eigen::Array3d position);
  double __interpolateVelocityV(Eigen::Array3d position);
  double __interpolateVelocityW(Eigen::Array3d position);
  double __interpolateVelocityU(Eigen::Array3d position, double &min,
                                double &max);
  double __interpolateVelocityV(Eigen::Array3d position, double &min,
                                double &max);
  double __interpolateVelocityW(Eigen::Array3d position, double &min,
                                double &max);
  double __interpolatePhi(Eigen::Array3d position, double &_min, double &_max);
  double __interpolatePhi(Eigen::Array3d position);

  virtual void print() {
    auto &phi = getScalarData(_phiId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/redistance/3d/" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);

    for (size_t i = 0; i < cellCount(); i++) {
      file << phi[i] << "\n";
    }
    file.close();
  }

protected:
  size_t _velocityId, _phiId;
  bool _isPressure2nd;
  double _dt, _tolerance;
};

} // namespace Leviathan

#endif