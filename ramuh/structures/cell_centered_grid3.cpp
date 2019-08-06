#include <structures/cell_centered_grid3.h>

namespace Ramuh {
CellCenteredGrid3::CellCenteredGrid3()
    : CellCenteredGrid3(
          BoundingBox3(Eigen::Array3d(-1, -1, -1), Eigen::Array3d(1, 1, 1)),
          Eigen::Array3i(32, 32, 32)) {}

CellCenteredGrid3::CellCenteredGrid3(BoundingBox3 domain,
                                     Eigen::Array3i gridSize)
    : _domain(domain), _gridSize(gridSize) {}

size_t CellCenteredGrid3::ijkToid(size_t i, size_t j, size_t k) {
  return k * _gridSize[0] * _gridSize[1] + j * _gridSize[0] + i;
}

std::tuple<size_t, size_t, size_t> CellCenteredGrid3::idToijk(size_t id) {
  size_t first, second, third;
  third = id / (_gridSize[0] * _gridSize[1]);
  second = (id % (_gridSize[0] * _gridSize[1])) / _gridSize[0];
  first = (id % (_gridSize[0] * _gridSize[1])) % _gridSize[0];
  return std::make_tuple(first, second, third);
}

size_t CellCenteredGrid3::newLabel(std::string label) {
  return newLabel(label, 0);
}

size_t CellCenteredGrid3::newLabel(std::string label, double value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _data.emplace_back(std::vector<double>(_gridSize.prod(), value));
    _dataLabel[label] = _data.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<double> &CellCenteredGrid3::getLabel(std::string label) {
  return _data[_dataLabel[label]];
}

double CellCenteredGrid3::operator()(std::string label, size_t i, size_t j,
                                     size_t k) {
  return _data[_dataLabel[label]][ijkToid(i, j, k)];
}

double CellCenteredGrid3::operator()(std::string label, size_t id) {
  return _data[_dataLabel[label]][id];
}

double CellCenteredGrid3::operator()(size_t i, size_t j, size_t k) {
  return _data[0][ijkToid(i, j, k)];
}

double CellCenteredGrid3::operator()(size_t id) { return _data[0][id]; }

} // namespace Ramuh