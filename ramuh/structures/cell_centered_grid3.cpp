#include <structures/cell_centered_grid3.h>
#include <iostream>

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

Eigen::Array3d CellCenteredGrid3::getH() {
  return _domain.size().cwiseQuotient(_gridSize.cast<double>());
}

Eigen::Array3d CellCenteredGrid3::getPosition(int i, int j, int k) {
  auto h = getH();
  return _domain.min() +
         Eigen::Array3d((i + 0.5) * h[0], (j + 0.5) * h[1], (k + 0.5) * h[2]);
}

Eigen::Array3d CellCenteredGrid3::getPosition(int id) {
  auto ijk = idToijk(id);
  return getPosition(ijk[0], ijk[1], ijk[2]);
}

std::vector<int> CellCenteredGrid3::idToijk(size_t id) {
  std::vector<int> indices(3);
  indices[2] = id / (_gridSize[0] * _gridSize[1]);
  indices[1] = (id % (_gridSize[0] * _gridSize[1])) / _gridSize[0];
  indices[0] = (id % (_gridSize[0] * _gridSize[1])) % _gridSize[0];
  return indices;
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

std::vector<double> &CellCenteredGrid3::getScalarData(std::string label) {
  if (_dataLabel.find(label) != _dataLabel.end())
    return _scalarData[_dataLabel[label]];
  std::cerr << "\033[1;31m[ ERROR ]\033[0mLabel not found!\n";
  return _scalarData[0];
}

std::vector<double> &CellCenteredGrid3::getScalarData(size_t index) {
  return _scalarData[index];
}

std::vector<Eigen::Array3d> &
CellCenteredGrid3::getArrayData(std::string label) {
  return _arrayData[_dataLabel[label]];
}

std::vector<Eigen::Array3d> &CellCenteredGrid3::getArrayData(size_t index) {
  return _arrayData[index];
}

size_t CellCenteredGrid3::cellCount() { return _gridSize.prod(); }

} // namespace Ramuh