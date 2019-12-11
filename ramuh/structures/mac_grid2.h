#ifndef __RAMUH_MAC_GRID_2_H__
#define __RAMUH_MAC_GRID_2_H__

#include <structures/cell_centered_grid2.h>
#include <geometry/vector2.h>

namespace Ramuh {
class MacGrid2 : public CellCenteredGrid2 {
public:
  MacGrid2();

  MacGrid2(BoundingBox2 domain, Eigen::Array2i gridSize);

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
   * @brief Returns a face position over the domain. The face coordinate should
   * be given first and then the face id which position is queried
   *
   * @param face x, y or z coordinate (0, 1 or 2)
   * @param faceId
   * @return Eigen::Array2d position in domain
   */
  Eigen::Array2d getFacePosition(size_t face, int faceId);
  Eigen::Array2d getFacePosition(size_t face, int i, int j);

  /**
   * @brief Return the total amount of faces over a single dimension, based on
   * the number of cells and its dimensions
   *
   * @return int total number of faces
   */
  int faceCount(int face);

  int faceijToid(int face, int i, int j);
  std::vector<int> faceIdToij(int face, int id);

  /**
   * @brief Get a reference for the vector containing data for a given face.
   *Face id should follow (x = 0, y = 1)
   *
   * @param face coordinate
   * @param which data are being retrieved
   * @return std::vector<double>& reference pointer to the data
   **/
  std::vector<double> &getFaceScalarData(size_t face, std::string label);
  std::vector<double> &getFaceScalarData(size_t face, int id);
  std::vector<Eigen::Array2d> &getFaceArrayData(size_t face, std::string label);
  std::vector<Eigen::Array2d> &getFaceArrayData(size_t face, int id);

protected:
  std::vector<std::vector<double>> _uScalar, _vScalar;
  std::vector<std::vector<Eigen::Array2d>> _uArray, _vArray;
  std::map<std::string, size_t> _faceDataLabel;
};

} // namespace Ramuh

#endif