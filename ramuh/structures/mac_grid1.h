#ifndef __RAMUH_MAC_GRID_1_H__
#define __RAMUH_MAC_GRID_1_H__

#include <structures/cell_centered_grid1.h>

namespace Ramuh {
class MacGrid1 : public CellCenteredGrid1 {
public:
  MacGrid1();

  MacGrid1(BoundingBox1 domain, int gridSize);

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
  std::vector<double> &getFaceLabelData(std::string label);
  std::vector<double> &getFaceLabelData(int id);

protected:
  std::vector<std::vector<double>> _uData;
  std::map<std::string, size_t> _faceLabel;
};

} // namespace Ramuh

#endif