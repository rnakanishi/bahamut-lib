#ifndef __RAMUH_CELL_CENTERED_GRID_3_H__
#define __RAMUH_CELL_CENTERED_GRID_3_H__

#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <geometry/bounding_box.h>
#include <Eigen/Dense>

namespace Ramuh {

class CellCenteredGrid3 {
public:
  CellCenteredGrid3();
  CellCenteredGrid3(BoundingBox3 domain, Eigen::Array3i gridSize);

  size_t ijkToid(size_t i, size_t j, size_t k);

  std::vector<int> idToijk(size_t id);

  /**
   * @brief Computes the cell size using domain information and number of grid
   *cells per dimension
   *
   * @return An Eigen array containing cell spacing for each dimension
   **/
  Eigen::Array3d getH();

  /**
   * @brief Get the Cell Array3d position object of a given cell id, or given
   * cell coordinates.
   *
   * @param id or i. If only one parameter is given, then the coordinates are
   * computed automatically
   * @param j second index coordinate
   * @param k third index coordinate
   * @return Array3d cell id's center position
   */
  Eigen::Array3d getCellPosition(int id);
  Eigen::Array3d getCellPosition(int i, int j, int k);

  /**
   * @brief Get the Cell Bounding Box object of a given cell id, or given cell
   * coordinates.
   *
   * @param id or i. If only one parameter is given, then the coordinates are
   * computed automatically
   * @param j second index coordinate
   * @param k third index coordinate
   * @return BoundingBox3 bounding box of the queried cell
   */
  BoundingBox3 getCellBoundingBox(int id);
  BoundingBox3 getCellBoundingBox(int i, int j, int k);

  /**
   * @brief Create a new scalar label in the structure. A initial value for the
   *grid can be assigned to all cells as well. If label already exists, then the
   *internal id of that label is returned.
   *
   * @param label for semantic purposes. Labels have to be unique over the same
   *instance
   * @param [optional] intialValue of the vecttor for that label
   * @return size_t return internal position of the label.
   **/
  size_t newCellScalarLabel(std::string label);
  size_t newCellScalarLabel(std::string label, double initialValue);

  size_t newCellArrayLabel(std::string label);
  size_t newCellArrayLabel(std::string label, Eigen::Array3d initialValue);

  /**
   * @brief Get the scalar vector object for a given label. This vector contains
   *whole data for that label.
   *
   * @param label string value. Must have been created
   * @return std::vector<double>& vector containing the data for that label
   **/
  std::vector<double> &getCellScalarData(std::string label);
  std::vector<double> &getCellScalarData(size_t index);

  std::vector<Eigen::Array3d> &getCellArrayData(std::string label);
  std::vector<Eigen::Array3d> &getCellArrayData(size_t index);

  /**
   * @brief Return the total number of cells present in the grid
   *
   * @return total number of cells
   **/
  size_t cellCount();

  /**
   * @brief Given a cell field id and a point inside the domain, computes an
   * inteprolated value at that position using trilinear method. If the point is
   * outside the domain, an extrapolated (and not accurate) value is returned
   *
   * @param dataId cell field which value is wanted
   * @param position point in the space as target
   * @return double interpolated value at position
   */
  double interpolateCellScalarData(int dataId, Eigen::Array3d position);

  /**
   * @brief Given a cell field id and a point inside the domain, computes an
   * inteprolated array value at that position using trilinear method. If the
   * point is outside the domain, an extrapolated (and not accurate) value is
   * returned
   *
   * @param dataId cell field which value is wanted
   * @param position point in the space as target
   * @return double interpolated value at position
   */
  Eigen::Array3d interpolateCellArrayData(int dataId, Eigen::Array3d position);

protected:
  Eigen::Array3i _gridSize;
  BoundingBox3 _domain;

  std::vector<std::vector<double>> _scalarData;
  std::vector<std::vector<Eigen::Array3d>> _arrayData;
  std::map<std::string, size_t> _dataLabel;
};

} // namespace Ramuh

#endif
