#ifndef __RAMUH_GRID3_H__
#define __RAMUH_GRID3_H__

#include <geometry/vector3.h>
#include <geometry/vector2.h>
#include <utils/material.h>
#include <vector>
#include <Eigen/Dense>

namespace Ramuh {

class RegularGrid3 {
public:
  RegularGrid3();

  ///
  /// Grid size
  /// \return Vector3d
  Vector3d gridSize();

  ///
  /// Grid domain size
  /// \return Vector3d
  Vector3d domainSize();

  ///
  /// Grid integer resolution
  /// \return Vector3i
  Vector3i resolution();

  ///
  /// Return the grid spacing values
  /// \return Vector3d spacing
  Vector3d h();

  ///
  /// Return the total amount of cells present on the grid
  /// \return int
  int cellCount();

  int ijkToId(int i, int j, int k);

  Vector3i idToijk(int id);

  ///
  /// Change grid size to the newSize. This method also updates grid spacing
  /// values
  /// \param newSize Vector3d containing the new size values
  virtual void setSize(Vector3d newSize);

  ///
  /// Change grid resolution. This method also updates grid spacing values
  /// \param newResolution
  virtual void setResolution(Vector3i newResolution);

  ///
  /// Computes velocities divergent in cell center and then solves pressure
  /// Poisson equation. After, updates the velocity values.
  virtual void solvePressure() ;

  ///
  /// Check boundary velocities and set them to the solid velocity (free slip)
  void boundaryVelocities();

  ///
  /// Add gravity value (\f$-9.81 \frac{m}{s^2}\f$) to all vertical velocities
  void addGravity();

  ///
  /// Make grid velocity advection using a semi lagrangian mehtod
  void advectGridVelocity();

  void macComarckVelocityAdvection();

  ///
  /// Initialize velocities
  void setVelocity();

  void extrapolateVelocity();

  void printFaceVelocity();

protected:
  double _interpolateVelocityU(Eigen::Array3d position);
  double _interpolateVelocityV(Eigen::Array3d position);
  double _interpolateVelocityW(Eigen::Array3d position);

  ///
  /// Change values for spacing. This method is called only through
  /// setResolution(...) or setSize(...) methods
  /// \param newH
  void setH(Vector3d newH);

  double _dt;                      // Time step
  Vector3i _resolution;            // Number of cells in each dimension
  Vector3d _domainSize;            // domain size in units
  Vector3d _h;                     // Spacing between cells
  std::vector<char> _materialMask; // Either cell is fluid, air, solid

  // TODO: Change to Matrix2 type
  // TODO: use two matrices to improve performance
  std::vector<std::vector<std::vector<Vector3d>>> _u, _v,
      _w; // Velocity components stored on faces
  std::vector<std::vector<std::vector<Material::FluidMaterial>>>
      _material; // Wheter the cell is a fluid, solid or air
};

} // namespace Ramuh

#endif