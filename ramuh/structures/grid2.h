#ifndef __RAMUH_GRID2_H__
#define __RAMUH_GRID2_H__

#include <geometry/vector3.h>
#include <geometry/vector2.h>
#include <utils/material.h>
#include <vector>

namespace Ramuh {

class RegularGrid2 {
public:
  RegularGrid2();

  ///
  /// Grid size
  /// \return Vector2d
  Vector2d gridSize();

  ///
  /// Grid domain size
  /// \return Vector2d
  Vector2d domainSize();

  ///
  /// Grid integer resolution
  /// \return Vector2i
  Vector2i resolution();

  ///
  /// Return the grid spacing values
  /// \return Vector2d spacing
  Vector2d h();

  ///
  /// Return the total amount of cells present on the grid
  /// \return int
  int cellCount();

  int ijToId(int i, int j);

  Vector2i idToij(int id);

  ///
  /// Change grid size to the newSize. This method also updates grid spacing
  /// values
  /// \param newSize Vector2d containing the new size values
  virtual void setSize(Vector2d newSize);

  ///
  /// Change grid resolution. This method also updates grid spacing values
  /// \param newResolution
  virtual void setResolution(Vector2i newResolution);

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

  ///
  /// Initialize velocities
  void setVelocity();

  void extrapolateVelocity();

  void printFaceVelocity();

protected:
  ///
  /// Change values for spacing. This method is called only through
  /// setResolution(...) or setSize(...) methods
  /// \param newH
  void setH(Vector2d newH);

  double _dt;                      // Time step
  Vector2i _resolution;            // Number of cells in each dimension
  Vector2d _domainSize;            // domain size in units
  Vector2d _h;                     // Spacing between cells
  std::vector<char> _materialMask; // Either cell is fluid, air, solid

  // TODO: Change to Matrix2 type
  // TODO: use two matrices to improve performance
  std::vector<std::vector<Vector2d>> _u,
      _v; // Velocity components stored on faces
  std::vector<std::vector<Material::FluidMaterial>>
      _material; // Wheter the cell is a fluid, solid or air
};

} // namespace Ramuh

#endif