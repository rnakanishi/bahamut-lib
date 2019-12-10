#include <structures/cell_centered_grid2.h>
#include <blas/interpolator.h>

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

Eigen::Array2i CellCenteredGrid2::getGridSize() { return _gridSize; }

Eigen::Array2i CellCenteredGrid2::getResolution() { return _gridSize; }

BoundingBox2 CellCenteredGrid2::getDomain() { return _domain; }

Eigen::Array2d CellCenteredGrid2::getH() {
  return _domain.getSize().cwiseQuotient(_gridSize.cast<double>());
}

Eigen::Array2d CellCenteredGrid2::getCellPosition(int i, int j) {
  auto h = getH();
  return _domain.getMin() + Eigen::Array2d((i + 0.5) * h[0], (j + 0.5) * h[1]);
}

Eigen::Array2d CellCenteredGrid2::getCellPosition(int id) {
  auto ij = idToij(id);
  auto h = getH();
  int i = ij[0], j = ij[1];
  return _domain.getMin() + Eigen::Array2d((i + 0.5) * h[0], (j + 0.5) * h[1]);
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

double CellCenteredGrid2::interpolateCellScalarData(int dataId,
                                                    Eigen::Array2d position) {
  auto &data = getCellScalarData(dataId);
  auto h = getH();

  // Find which cell-id data belongs
  Eigen::Array2i cellId =
      (position - _domain.getMin()).cwiseQuotient(h).floor().cast<int>();

  // Assemble bilinear stencil interpolation for velocities
  auto cellPos = getCellPosition(cellId[0], cellId[1]);

  int xindex, yindex;
  xindex = std::min(_gridSize[0] - 1, cellId[0] + 1);
  yindex = std::min(_gridSize[1] - 1, cellId[1] + 1);
  if (position[0] < cellPos[0]) {
    xindex = std::max(0, cellId[0] - 1);
    h[0] = -h[0];
  }
  if (position[1] < cellPos[1]) {
    yindex = std::max(0, cellId[1] - 1);
    h[1] = -h[1];
  }

  std::vector<Eigen::Array2d> points;
  std::vector<double> values(4);
  double target[2];
  target[0] = position[0];
  target[1] = position[1];

  // cell Stencil for linear interpolation
  points.emplace_back(cellPos);
  points.emplace_back(cellPos + Eigen::Array2d(h[0], 0));
  points.emplace_back(cellPos + Eigen::Array2d(0, h[1]));
  points.emplace_back(cellPos + h);

  values[0] = data[ijToid(cellId[0], cellId[1])];
  values[1] = data[ijToid(xindex, cellId[1])];
  values[2] = data[ijToid(cellId[0], yindex)];
  values[3] = data[ijToid(xindex, yindex)];

  return Ramuh::Interpolator::bilinear(target, points, values);
}

Eigen::Array2d
CellCenteredGrid2::interpolateCellArrayData(int dataId,
                                            Eigen::Array2d position) {
  auto &data = getCellArrayData(dataId);
  auto h = getH();

  // Find which cell-id data belongs
  Eigen::Array2i cellId =
      (position - _domain.getMin()).cwiseQuotient(h).floor().cast<int>();

  // Assemble bilinear stencil interpolation for velocities
  auto cellPos = getCellPosition(cellId[0], cellId[1]);

  int xindex, yindex;
  xindex = std::min(_gridSize[0] - 1, cellId[0] + 1);
  yindex = std::min(_gridSize[1] - 1, cellId[1] + 1);
  if (position[0] < cellPos[0]) {
    xindex = std::max(0, cellId[0] - 1);
    h[0] = -h[0];
  }
  if (position[1] < cellPos[1]) {
    yindex = std::max(0, cellId[1] - 1);
    h[1] = -h[1];
  }

  std::vector<Eigen::Array2d> points;
  std::vector<double> values(4);
  double target[2];
  target[0] = position[0];
  target[1] = position[1];

  Eigen::Array2d interpData;
  // cell Stencil for linear interpolation
  points.emplace_back(cellPos);
  points.emplace_back(cellPos + Eigen::Array2d(h[0], 0));
  points.emplace_back(cellPos + Eigen::Array2d(0, h[1]));
  points.emplace_back(cellPos + h);

  for (size_t d = 0; d < 2; d++) {
    values[0] = data[ijToid(cellId[0], cellId[1])][d];
    values[1] = data[ijToid(xindex, cellId[1])][d];
    values[2] = data[ijToid(cellId[0], yindex)][d];
    values[3] = data[ijToid(xindex, yindex)][d];
    interpData[d] = Ramuh::Interpolator::bilinear(target, points, values);
  }
  return interpData;
}

} // namespace Ramuh