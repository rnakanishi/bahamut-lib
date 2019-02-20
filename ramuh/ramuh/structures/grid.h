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
  /// Grid size
  /// \return Vector3d
  Vector3d size();

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
  void setSize(Vector3d newSize);

  ///
  /// Change grid resolution. This method also updates grid spacing values
  /// \param newResolution
  void setResolution(Vector3i newResolution);

protected:
  ///
  /// Change values for spacing. This method is called only through
  /// setResolution(...) or setSize(...) methods
  /// \param newH
  void setH(Vector3d newH);

  Vector3i _resolution;            // Number of cells in each dimension
  Vector3d _size;                  // domain size in units
  Vector3d _h;                     // Spacing between cells
  std::vector<char> _materialMask; // Either cell is fluid, air, solid
  std::vector<double> u, v, w;     // Velocity components stored on faces
};

} // namespace Ramuh