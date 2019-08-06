#ifndef __RAMUH_GRID2_H__
#define __RAMUH_GRID2_H__

#include <geometry/vector3.h>
#include <geometry/vector2.h>
#include <utils/material.h>
#include <vector>
#include <Eigen/Dense>

namespace Ramuh {

class RegularGrid2 {
public:
  RegularGrid2();

  ///
  /// Grid size
  /// \return Eigen::Array2d
  Eigen::Array2d gridSize();

  ///
  /// Grid domain size
  /// \return Eigen::Array2d
  Eigen::Array2d domainSize();

  ///
  /// Grid integer resolution
  /// \return Eigen::Array2i
  Eigen::Array2i resolution();

  ///
  /// Return the grid spacing values
  /// \return Eigen::Array2d spacing
  Eigen::Array2d h();

  ///
  /// Return the total amount of cells present on the grid
  /// \return int
  int cellCount();

  int ijToId(int i, int j);

  Eigen::Array2i idToij(int id);

  void setDt(double dt);

  ///
  /// Change grid size to the newSize. This method also updates grid spacing
  /// values
  /// \param newSize Eigen::Array2d containing the new size values
  virtual void setSize(Eigen::Array2d newSize);

  ///
  /// Change grid resolution. This method also updates grid spacing values
  /// \param newResolution
  virtual void setResolution(Eigen::Array2i newResolution);

  ///
  /// Computes velocities divergent in cell center and then solves pressure
  /// Poisson equation. After, updates the velocity values.
  virtual void solvePressure();

  ///
  /// Check boundary velocities and set them to the solid velocity (free slip)
  void boundaryVelocities();

  ///
  /// Add gravity value (\f$-9.81 \frac{m}{s^2}\f$) to all vertical velocities
  void addGravity();

  ///
  /// Make grid velocity advection using a semi lagrangian mehtod
  void advectGridVelocity();

  void macCormackVelocityAdvection();

  ///
  /// Initialize velocities
  void setVelocity();

  void extrapolateVelocity();

  void printFaceVelocity();

  void cfl();

  bool advanceTime();

protected:
  ;
  ;
  ;
  double _interpolateVelocityU(Eigen::Array2d position, double &min,
                               double &max);
  double _interpolateVelocityV(Eigen::Array2d position, double &min,
                               double &max);
  double _interpolateVelocityU(Eigen::Array2d position);
  double _interpolateVelocityV(Eigen::Array2d position);

  /**
   * @brief Return if both sides of a given cell have the same material
   *
   * @param cellId query cell
   * @param material Material we are querying
   * @return true if both cells have same material
   * @return false otherwise
   */
  bool _hasOppositeNeighborsWithMaterial(int cellId,
                                         Material::FluidMaterial material);
  ///
  /// Change values for spacing. This method is called only through
  /// setResolution(...) or setSize(...) methods
  /// \param newH
  void setH(Eigen::Array2d newH);

  double _tolerance;
  double _ellapsedDt;
  double _originalDt;
  double _dt;                      // Time step
  Eigen::Array2i _resolution;      // Number of cells in each dimension
  Eigen::Array2d _domainSize;      // domain size in units
  Eigen::Array2d _h;               // Spacing between cells
  std::vector<char> _materialMask; // Either cell is fluid, air, solid

  // TODO: Change to Matrix2 type
  // TODO: use two matrices to improve performance
  std::vector<std::vector<double>> _u,
      _v; // Velocity components stored on faces
  std::vector<std::vector<Material::FluidMaterial>>
      _material; // Wheter the cell is a fluid, solid or air
  std::vector<std::vector<Material::FluidMaterial>> _uFaceMaterial;
  std::vector<std::vector<Material::FluidMaterial>> _vFaceMaterial;
  std::vector<int> _fluidCells;
};

} // namespace Ramuh

#endif