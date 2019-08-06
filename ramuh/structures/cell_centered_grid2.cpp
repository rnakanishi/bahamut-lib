#include <structures/cell_centered_grid2.h>

namespace Ramuh {
CellCenteredGrid2::CellCenteredGrid2() : _gridSize(Eigen::Array2i(32, 32)) {}

void CellCenteredGrid2::setGridSize(Eigen::Array2i size) { _gridSize = size; }

void CellCenteredGrid2::setGridSize(std::pair<size_t, size_t> size) {
  setGridSize(Eigen::Array2i(size.first, size.second));
}

size_t CellCenteredGrid2::ijToid(size_t i, size_t j) {
  return j * _gridSize[0] + i;
}

std::pair<size_t, size_t> CellCenteredGrid2::idToij(size_t id) {
  std::pair<size_t, size_t> index;
  index.first = id / (_gridSize[0]);
  index.second = id % _gridSize[0];
  return index;
}

size_t CellCenteredGrid2::newLabel(std::string label) {
  return newLabel(label, 0);
}

size_t CellCenteredGrid2::newLabel(std::string label, double value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _data.emplace_back(std::vector<double>(_gridSize.prod(), value));
    _dataLabel[label] = _data.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<double> &CellCenteredGrid2::getLabel(std::string label) {
  return _data[_dataLabel[label]];
}

double CellCenteredGrid2::operator()(std::string label, size_t i, size_t j) {
  return _data[_dataLabel[label]][ijToid(i, j)];
}

double CellCenteredGrid2::operator()(std::string label, size_t id) {
  return _data[_dataLabel[label]][id];
}

double CellCenteredGrid2::operator()(size_t i, size_t j) {
  return _data[0][ijToid(i, j)];
}

double CellCenteredGrid2::operator()(size_t id) { return _data[0][id]; }

} // namespace Ramuh