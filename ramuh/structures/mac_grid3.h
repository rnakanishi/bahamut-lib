#ifndef __RAMUH_MAC_GRID_3_H__
#define __RAMUH_MAC_GRID_3_H__
#include <structures/cell_centered_grid3.h>

namespace Ramuh {
class MacGrid3 : public CellCenteredGrid3 {
public:
  MacGrid3();
  MacGrid3(BoundingBox3 domain, Eigen::Array3i gridSize);

  /**
   * @brief Creates a new label for the data to be stored. Labels should be
   *unique among the same grid structure. The methods return an id for the
   *created vector
   *
   * @param semantic label for the new data
   * @param [opt] initial value for all faces. Defalut is zero.
   * @return size_t
   **/
  size_t newFaceScalarLabel(std::string label);
  size_t newFaceScalarLabel(std::string label, double value);

  size_t newFaceArrayLabel(std::string label);
  size_t newFaceArrayLabel(std::string label, Eigen::Array3d value);

  int getFaceScalarLabelId(std::string label);

  int getFaceArrayLabelId(std::string label);

  /**
   * @brief Get a reference for the vector containing data for a given face.
   *Face id should follow (x = 0, y = 1, z = 2)
   *
   * @param face coordinate
   * @param which data are being retrieved
   * @return std::vector<double>& reference pointer to the data
   **/
  std::vector<double> &getFaceScalarData(size_t face, std::string label);
  std::vector<double> &getFaceScalarData(size_t face, size_t index);
  std::vector<Eigen::Array3d> &getFaceArrayData(size_t face, std::string label);
  std::vector<Eigen::Array3d> &getFaceArrayData(size_t face, size_t index);

  Eigen::Array3d getFacePosition(int face, int id);
  Eigen::Array3d getFacePosition(int face, int i, int j, int k);

  /**
   * @brief Face grid size differ from cell grid size due to values being placed
   * at faces. Each coordinate has one more value then the cell counterpart.
   *
   * @param face x, y or z face (0,1 or 2 respectively)
   * @return Eigen::Array3i Increments one into the grid size of corresponding
   * face
   */
  Eigen::Array3i getFaceGridSize(int face);

  int faceCount(int face);

  int faceijkToid(int face, int i, int j, int k);
  std::vector<int> faceIdToijk(int face, int id);

  /**
   * @brief
   *
   * @param dataId
   * @param position
   * @return double
   */
  double interpolateFaceScalarData(int face, int dataId,
                                   Eigen::Array3d position);

  /**
   * @brief
   *
   * @param dataId
   * @param position
   * @return Eigen::Array3d
   */
  Eigen::Array3d interpolateFaceArrayData(int face, int dataId,
                                          Eigen::Array3d position);

protected:
  std::vector<std::vector<double>> _uScalar, _vScalar, _wScalar;
  std::vector<std::vector<Eigen::Array3d>> _uArray, _vArray, _wArray;
  std::map<std::string, size_t> _faceScalarLabel, _faceArrayLabel;
};

} // namespace Ramuh

#endif