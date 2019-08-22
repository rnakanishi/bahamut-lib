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

std::pair<size_t, size_t> CellCenteredGrid2::idToij(size_t id) {
  std::pair<size_t, size_t> index;
  index.second = id / (_gridSize[0]);
  index.first = id % _gridSize[0];
  return index;
}

Eigen::Array2d CellCenteredGrid2::getH() {
  return _domain.size().cwiseQuotient(_gridSize.cast<double>());
}

Eigen::Array2d CellCenteredGrid2::getPosition(int i, int j) {
  auto h = getH();
  return _domain.min() + Eigen::Array2d((i + 0.5) * h[0], (j + 0.5) * h[1]);
}

Eigen::Array2d CellCenteredGrid2::getPosition(int id) {
  auto ij = idToij(id);
  auto h = getH();
  int i = ij.first, j = ij.second;
  return _domain.min() + Eigen::Array2d((i + 0.5) * h[0], (j + 0.5) * h[1]);
}

int CellCenteredGrid2::cellCount() { return _gridSize.prod(); }

size_t CellCenteredGrid2::newScalarLabel(std::string label) {
  return newScalarLabel(label, 0);
}

size_t CellCenteredGrid2::newScalarLabel(std::string label, double value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _scalarData.emplace_back(std::vector<double>(_gridSize.prod(), value));
    _dataLabel[label] = _scalarData.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<double> &CellCenteredGrid2::getScalarLabel(std::string label) {
  return _scalarData[_dataLabel[label]];
}

std::vector<double> &CellCenteredGrid2::getScalarLabel(size_t id) {
  return _scalarData[id];
}

size_t CellCenteredGrid2::newArrayLabel(std::string label) {
  return newArrayLabel(label, Eigen::Array2d(0., 0.));
}

size_t CellCenteredGrid2::newArrayLabel(std::string label,
                                        Eigen::Array2d value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _arrayData.emplace_back(
        std::vector<Eigen::Array2d>(_gridSize.prod(), value));
    _dataLabel[label] = _arrayData.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<Eigen::Array2d> &
CellCenteredGrid2::getArrayLabel(std::string label) {
  return _arrayData[_dataLabel[label]];
}

} // namespace Ramuh