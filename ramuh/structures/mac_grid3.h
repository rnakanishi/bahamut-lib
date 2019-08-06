#ifndef __RAMUH_MAC_GRID_H__
#define __RAMUH_MAC_GRID_H__

#include <structures/cell_centered_grid3.h>

namespace Ramuh {
class MacGrid3 : CellCenteredGrid3 {
public:
  MacGrid3();

  /**
   * @brief Creates a new label for the data to be stored. Labels should be
   *unique among the same grid structure. The methods return an id for the
   *created vector
   *
   * @param semantic label for the new data
   * @param [opt] initial value for all faces. Defalut is zero.
   * @return size_t
   **/
  size_t newFaceLabel(std::string label);
  size_t newFaceLabel(std::string label, double value);

  /**
   * @brief Get a reference for the vector containing data for a given face.
   *Face id should follow (x = 0, y = 1, z = 2)
   *
   * @param face coordinate
   * @param which data are being retrieved
   * @return std::vector<double>& reference pointer to the data
   **/
  std::vector<double> &getFaceLabel(size_t face, std::string label);

private:
  std::vector<std::vector<double>> _udata, _vdata, _wdata;
  std::map<std::string, size_t> _faceDataLabel;
};

} // namespace Ramuh

#endif