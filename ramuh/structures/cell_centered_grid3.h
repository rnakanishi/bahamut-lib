#ifndef __RAMUH_CELL_CENTERED_GRID_H__
#define __RAMUH_CELL_CENTERED_GRID_H__

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

  std::tuple<size_t, size_t, size_t> idToijk(size_t id);

  /**
   * @brief Create a new label in the structure. A initial value for the grid
   *can be assigned to all cells as well. If label already exists, then the
   *internal id of that label is returned.
   *
   * @param label for semantic purposes. Labels have to be unique over the same
   *instance
   * @return size_t return internal position of the label.
   **/
  size_t newLabel(std::string label);
  size_t newLabel(std::string label, double initialValue);

  /**
   * @brief Get the vector object for a given label. This vector contains whole
   *data for that label.
   *
   * @param label string value. Must have been created
   * @return std::vector<double>& vector containing the data for that label
   **/
  std::vector<double> &getLabel(std::string label);

  /**
   * @brief This operator access the data from a given label and given
   *coordinates, or cell id. If no label is assigned, the the first created
   *label is accessed
   *
   * @param label label to access
   * @param if only one value is assigned, then call using the cell id value.
   *Otherwise, calls using the ijk coordinate
   * @return double value of the corresponding label-coordinate
   **/
  double operator()(std::string label, size_t i, size_t j, size_t k);
  double operator()(std::string label, size_t id);
  double operator()(size_t i, size_t j, size_t k);
  double operator()(size_t id);

protected:
  Eigen::Array3i _gridSize;
  BoundingBox3 _domain;

  std::vector<std::vector<double>> _data;
  std::map<std::string, size_t> _dataLabel;
};

} // namespace Ramuh

#endif
