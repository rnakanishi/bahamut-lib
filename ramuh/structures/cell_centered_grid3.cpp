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

size_t CellCenteredGrid3::newScalarLabel(std::string label) {
  return newScalarLabel(label, 0);
}

size_t CellCenteredGrid3::newScalarLabel(std::string label, double value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _scalarData.emplace_back(std::vector<double>(_gridSize.prod(), value));
    _dataLabel[label] = _scalarData.size() - 1;
  }
  return _dataLabel[label];
}

size_t CellCenteredGrid3::newArrayLabel(std::string label) {
  return newArrayLabel(label, Eigen::Array3d(0, 0, 0));
}

size_t CellCenteredGrid3::newArrayLabel(std::string label,
                                        Eigen::Array3d value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _arrayData.emplace_back(
        std::vector<Eigen::Array3d>(_gridSize.prod(), value));
    _dataLabel[label] = _arrayData.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<double> &CellCenteredGrid3::getScalarLabel(std::string label) {
  return _scalarData[_dataLabel[label]];
}

std::vector<Eigen::Array3d> &
CellCenteredGrid3::getArrayLabel(std::string label) {
  return _arrayData[_dataLabel[label]];
}

} // namespace Ramuh