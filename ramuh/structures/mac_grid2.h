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
  size_t newFaceScalarLabel(std::string label);
  size_t newFaceScalarLabel(std::string label, double value);

  size_t newFaceArrayLabel(std::string label);
  size_t newFaceArrayLabel(std::string label, Eigen::Array2d value);

  /**
   * @brief Get a reference for the vector containing data for a given face.
   *Face id should follow (x = 0, y = 1)
   *
   * @param face coordinate
   * @param which data are being retrieved
   * @return std::vector<double>& reference pointer to the data
   **/
  std::vector<double> &getFaceScalarLabel(size_t face, std::string label);
  std::vector<Eigen::Array2d> &getFaceArrayLabel(size_t face,
                                                 std::string label);

protected:
  std::vector<std::vector<double>> _uScalar, _vScalar;
  std::vector<std::vector<Eigen::Array2d>> _uArray, _vArray;
  std::map<std::string, size_t> _faceDataLabel;
};

} // namespace Ramuh

#endif