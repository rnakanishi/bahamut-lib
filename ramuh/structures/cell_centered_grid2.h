#ifndef __RAMUH_CELL_CENTERED_GRID_2_H__
#define __RAMUH_CELL_CENTERED_GRID_2_H__

#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <geometry/bounding_box.h>
#include <Eigen/Dense>

namespace Ramuh {

class CellCenteredGrid2 {
public:
  CellCenteredGrid2();
  CellCenteredGrid2(BoundingBox2 domain, Eigen::Array2i gridSize);

  /**
   * @brief Set the Grid Size in number of cells in each coordinate (x,y)
   *
   * @param size new size to be used.
   **/
  void setGridSize(std::pair<size_t, size_t> size);
  void setGridSize(Eigen::Array2i size);

  /**
   * @brief From a pair of coordinate values (i,j), computes the single id value
   *using the grid size provided ealier
   *
   * @param i x coordinate value
   * @param j y coordinate value
   * @return size_t corresponding  id value
   **/
  size_t ijToid(size_t i, size_t j);

  /**
   * @brief Given an id value, computes the pair (i,j) of values that
   *generated it
   *
   * @param id
   * @return std::pair<size_t, size_t>
   **/
  std::vector<size_t> idToij(size_t id);

  Eigen::Array2d getH();

  Eigen::Array2d cellPosition(int i, int j);
  Eigen::Array2d cellPosition(int id);

  int cellCount();

  /**
   * @brief Create a new label in the structure. A initial value for the grid
   *can be assigned to all cells as well. If no initial value is given, then all
   *values are set to zero. If label already exists, then the internal id of
   *that label is returned.
   *
   * @param label for semantic purposes. Labels have to be unique over the same
   *instance
   * @return size_t return internal position of the label.
   **/
  size_t newCellScalarLabel(std::string label);
  size_t newCellScalarLabel(std::string label, double initialValue);

  size_t newCellArrayLabel(std::string label);
  size_t newCellArrayLabel(std::string label, Eigen::Array2d initialValue);

  /**
   * @brief Get the vector object for a given label. This vector contains whole
   *data for that label.
   *
   * @param label string value. Must have been created
   * @return std::vector<double>& vector containing the data for that label
   **/
  std::vector<double> &getCellScalarData(std::string label);
  std::vector<double> &getCellScalarData(size_t id);

  std::vector<Eigen::Array2d> &getCellArrayData(std::string label);
  std::vector<Eigen::Array2d> &getCellArrayData(int id);

protected:
  Eigen::Array2i _gridSize;
  BoundingBox2 _domain;

  std::vector<std::vector<double>> _scalarData;
  std::vector<std::vector<Eigen::Array2d>> _arrayData;
  std::map<std::string, size_t> _dataLabel;
};

} // namespace Ramuh

#endif
