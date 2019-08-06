#ifndef __RAMUH_MAC_GRID_H__
#define __RAMUH_MAC_GRID_H__

#include <structures/cell_centered_grid2.h>

namespace Ramuh {
class MacGrid2 : public CellCenteredGrid2 {
public:
  MacGrid2();

  /**
   * @brief Creates a new label for the data to be stored. Labels should be
   *unique among the same grid structure. The methods return an id for the
   *created vector
   *
   * @param label for the new semantic data
   * @param [opt] initial value for all faces
   * @return size_t
   **/
  size_t newFaceLabel(std::string label);
  size_t newFaceLabel(std::string label, double value);

  /**
   * @brief Get a reference for the vector containing data for a given face.
   *Face id should follow (x = 0, y = 1)
   *
   * @param face coordinate
   * @param which data are being retrieved
   * @return std::vector<double>& reference pointer to the data
   **/
  std::vector<double> &getFaceLabel(size_t face, std::string label);

private:
  std::vector<std::vector<double>> _udata, _vdata;
  std::map<std::string, size_t> _faceDataLabel;
};

} // namespace Ramuh

#endif