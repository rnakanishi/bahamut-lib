#ifndef __RAMUH_GRID_H__
#define __RAMUH_GRID_H__

#include <structures/vector3.h>
#include <vector>

namespace Ramuh {

class RegularGrid {
public:
  RegularGrid();

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
  /// Change grid size to the newSize. This method also updates grid spacing
  /// values
  /// \param newSize Vector3d containing the new size values
  virtual void setSize(Vector3d newSize);

  ///
  /// Change grid resolution. This method also updates grid spacing values
  /// \param newResolution
  virtual void setResolution(Vector3i newResolution);

  void setVelocity();

  void printFaceVelocity();

protected:
  ///
  /// Change values for spacing. This method is called only through
  /// setResolution(...) or setSize(...) methods
  /// \param newH
  void setH(Vector3d newH);

  Vector3i _resolution;            // Number of cells in each dimension
  Vector3d _domainSize;            // domain size in units
  Vector3d _h;                     // Spacing between cells
  std::vector<char> _materialMask; // Either cell is fluid, air, solid
  std::vector<std::vector<std::vector<double>>> _u, _v,
      _w; // Velocity components stored on faces
  std::vector<std::vector<std::vector<double>>>
      _pressure; // pressure component stored on cell center
};

} // namespace Ramuh

#endif