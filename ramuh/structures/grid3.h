#ifndef __RAMUH_GRID3_H__
#define __RAMUH_GRID3_H__

#include <geometry/vector3.h>
#include <geometry/vector2.h>
#include <utils/material.h>
#include <vector>
#include <Eigen/Dense>
#include <utility>

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

  double tolerance();

  void setTolerance(double tol);
  

  ///
  /// Return the total amount of cells present on the grid
  /// \return int
  int cellCount();

  /**
   * @brief Return faces that belong to the cell accordingly to the coordinate
   *
   * @param cell cell id
   * @param coord x y, or z coordinate
   * @return std::pair<Eigen::Array3i, Eigen::Array3i>
   **/
  std::vector<Eigen::Array3i> cellFaces(Eigen::Array3i cell, int coord);

  int ijkToId(int i, int j, int k);

  Eigen::Array3i idToijk(int id);

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
  virtual void solvePressure();

  ///
  /// Check boundary velocities and set them to the solid velocity (free slip)
  void boundaryVelocities();

  ///
  /// Add gravity value (\f$-9.81 \frac{m}{s^2}\f$) to all vertical velocities
  void addGravity();

  void addExternalForce(Eigen::Vector3d force);

  ///
  /// Make grid velocity advection using a semi lagrangian mehtod
  void advectGridVelocity();

  void macComarckVelocityAdvection();

  ///
  /// Initialize velocities
  void setVelocity();

  void extrapolateVelocity();

  void writeFaceVelocity(const char *filename);

  void readFaceVelocity(const char *filename);

  void cfl();

  bool advanceTime();

  void swapBuffers();

protected:
  double _interpolateVelocityU(Eigen::Array3d position);
  double _interpolateVelocityV(Eigen::Array3d position);
  double _interpolateVelocityW(Eigen::Array3d position);
  double _interpolateVelocityU(Eigen::Array3d position, double &_min,
                               double &_max);
  double _interpolateVelocityV(Eigen::Array3d position, double &_min,
                               double &_max);
  double _interpolateVelocityW(Eigen::Array3d position, double &_min,
                               double &_max);

  ///
  /// Change values for spacing. This method is called only through
  /// setResolution(...) or setSize(...) methods
  /// \param newH
  void setH(Vector3d newH);

  double _tolerance;
  double _maxVelocity[3];
  double _dt;           // CFL max timestep
  double _ellapsedDt;   // Total ellapsed time in frame
  double _originalDt;   // Time step
  Vector3i _resolution; // Number of cells in each dimension
  Vector3d _domainSize; // domain size in units
  Vector3d _h;          // Spacing between cells

  // TODO: Change to Matrix2 type
  // TODO: use two matrices to improve performance
  int _currBuffer;
  std::vector<std::vector<std::vector<double>>> _u[2], _v[2],
      _w[2]; // Velocity components stored on faces
  std::vector<int> _fluidCells, _surfaceCells;
  std::vector<std::vector<std::vector<Material::FluidMaterial>>>
      _material; // Wheter the cell is a fluid, solid or air
};

} // namespace Ramuh

#endif