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

  void advectEuler();

  void advectRungeKutta3();

  void advectMacCormack();

  void advectWeno();

  void advectCip();

  void advectUpwind();

  void computeCellsGradient();

  void computeWenoGradient();

  void computeCellVelocity();

  void redistance();

  bool advanceTime();

  void applyCfl();

  /**
   * @brief Find all cells that are close to the interface. This is made
   * computing the product between two cells. If the resulting value is
   * negative, than the cell is possibly an interface cell. All the cells are
   * internally marked and also their ids are returned in a std::vector
   * container
   *
   * @return std::vector<int> contains all the cell ids that are part of the
   * interface
   */
  std::vector<int> trackSurface();

protected:
  virtual void print() {
    auto &phi = getCellScalarData(_phiId);
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
  size_t _cellVelocityId, _phiId, _cellGradientId;
  size_t _faceVelocityId;
  bool _isPressure2nd;
  double _dt, _originalDt, _ellapsedDt;
  double _tolerance;

  /// All cells that are part of the interface are marked as true
  std::vector<bool> _surfaceCells;
};

} // namespace Leviathan

#endif