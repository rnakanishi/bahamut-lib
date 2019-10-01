#include <structures/cell_centered_grid2.h>

namespace Ramuh {
CellCenteredGrid2::CellCenteredGrid2()
    : CellCenteredGrid2(
          BoundingBox2(Eigen::Array2d(-1, -1), Eigen::Array2d(1, 1)),
          Eigen::Array2i(32, 32)) {}

CellCenteredGrid2::CellCenteredGrid2(BoundingBox2 domain,
                                     Eigen::Array2i gridSize)
    : _domain(domain), _gridSize(gridSize) {}

void CellCenteredGrid2::setGridSize(Eigen::Array2i size) { _gridSize = size; }

void CellCenteredGrid2::setGridSize(std::pair<size_t, size_t> size) {
  setGridSize(Eigen::Array2i(size.first, size.second));
}

size_t CellCenteredGrid2::ijToid(size_t i, size_t j) {
  return j * _gridSize[0] + i;
}

std::vector<size_t> CellCenteredGrid2::idToij(size_t id) {
  std::vector<size_t> index(2);
  index[0] = id % (_gridSize[0]);
  index[1] = id / _gridSize[0];
  return index;
}

Eigen::Array2d CellCenteredGrid2::getH() {
  return _domain.size().cwiseQuotient(_gridSize.cast<double>());
}

Eigen::Array2d CellCenteredGrid2::getCellPosition(int i, int j) {
  auto h = getH();
  return _domain.min() + Eigen::Array2d((i + 0.5) * h[0], (j + 0.5) * h[1]);
}

Eigen::Array2d CellCenteredGrid2::getCellPosition(int id) {
  auto ij = idToij(id);
  auto h = getH();
  int i = ij[0], j = ij[1];
  return _domain.min() + Eigen::Array2d((i + 0.5) * h[0], (j + 0.5) * h[1]);
}

BoundingBox2 CellCenteredGrid2::getCellBoundingBox(int i, int j) {
  auto position = getCellPosition(i, j);
  BoundingBox2 box;
  Eigen::Array2d h = getH();
  box.setMin(position - h / 2);
  box.setMax(position + h / 2);
  return box;
}

BoundingBox2 CellCenteredGrid2::getCellBoundingBox(int id) {
  auto ij = idToij(id);
  return getCellBoundingBox(ij[0], ij[1]);
}

int CellCenteredGrid2::cellCount() { return _gridSize.prod(); }

size_t CellCenteredGrid2::newCellScalarLabel(std::string label) {
  return newCellScalarLabel(label, 0);
}

size_t CellCenteredGrid2::newCellScalarLabel(std::string label, double value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _scalarData.emplace_back(std::vector<double>(_gridSize.prod(), value));
    _dataLabel[label] = _scalarData.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<double> &CellCenteredGrid2::getCellScalarData(std::string label) {
  return _scalarData[_dataLabel[label]];
}

std::vector<double> &CellCenteredGrid2::getCellScalarData(size_t id) {
  return _scalarData[id];
}

size_t CellCenteredGrid2::newCellArrayLabel(std::string label) {
  return newCellArrayLabel(label, Eigen::Array2d(0., 0.));
}

size_t CellCenteredGrid2::newCellArrayLabel(std::string label,
                                            Eigen::Array2d value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _arrayData.emplace_back(
        std::vector<Eigen::Array2d>(_gridSize.prod(), value));
    _dataLabel[label] = _arrayData.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<Eigen::Array2d> &
CellCenteredGrid2::getCellArrayData(std::string label) {
  return _arrayData[_dataLabel[label]];
}

std::vector<Eigen::Array2d> &CellCenteredGrid2::getCellArrayData(int id) {
  return _arrayData[id];
}

} // namespace Ramuh