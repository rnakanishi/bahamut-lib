#include <structures/cell_centered_grid1.h>

namespace Ramuh {
CellCenteredGrid1::CellCenteredGrid1()
    : CellCenteredGrid1(BoundingBox1(-1, 1), 32) {}

CellCenteredGrid1::CellCenteredGrid1(BoundingBox1 domain, int gridSize)
    : _domain(domain), _gridSize(gridSize) {}

void CellCenteredGrid1::setGridSize(int size) { _gridSize = size; }

double CellCenteredGrid1::getH() { return _domain.size() / _gridSize; }

double CellCenteredGrid1::getPosition(int i) {
  auto h = getH();
  return _domain.min() + (i + 0.5) * h;
}

size_t CellCenteredGrid1::newLabel(std::string label) {
  return newLabel(label, 0);
}

size_t CellCenteredGrid1::newLabel(std::string label, double value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _scalarData.emplace_back(std::vector<double>(_gridSize, value));
    _dataLabel[label] = _scalarData.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<double> &CellCenteredGrid1::getLabelData(std::string label) {
  return _scalarData[_dataLabel[label]];
}

std::vector<double> &CellCenteredGrid1::getLabelData(size_t id) {
  return _scalarData[id];
}

} // namespace Ramuh