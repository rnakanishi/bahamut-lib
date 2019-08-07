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
   * @brief Computes the cell size using domain information and number of grid
   *cells per dimension
   *
   * @return An Eigen array containing cell spacing for each dimension
   **/
  Eigen::Array3d getH();

  /**
   * @brief Create a new scalar label in the structure. A initial value for the
   *grid can be assigned to all cells as well. If label already exists, then the
   *internal id of that label is returned.
   *
   * @param label for semantic purposes. Labels have to be unique over the same
   *instance
   * @return size_t return internal position of the label.
   **/
  size_t newScalarLabel(std::string label);
  size_t newScalarLabel(std::string label, double initialValue);

  size_t newArrayLabel(std::string label);
  size_t newArrayLabel(std::string label, Eigen::Array3d initialValue);

  /**
   * @brief Get the scalar vector object for a given label. This vector contains
   *whole data for that label.
   *
   * @param label string value. Must have been created
   * @return std::vector<double>& vector containing the data for that label
   **/
  std::vector<double> &getScalarVector(std::string label);
  std::vector<double> &getScalarVector(size_t index);

  std::vector<Eigen::Array3d> &getArrayVector(std::string label);
  std::vector<Eigen::Array3d> &getArrayVector(size_t index);

protected:
  Eigen::Array3i _gridSize;
  BoundingBox3 _domain;

  std::vector<std::vector<double>> _scalarData;
  std::vector<std::vector<Eigen::Array3d>> _arrayData;
  std::map<std::string, size_t> _dataLabel;
};

} // namespace Ramuh

#endif
