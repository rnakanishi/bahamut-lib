#include <structures/cell_centered_grid3.h>
#include <blas/interpolator.h>

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
  return _domain.getSize().cwiseQuotient(_gridSize.cast<double>());
}

Eigen::Array3d CellCenteredGrid3::getCellPosition(int i, int j, int k) {
  auto h = getH();
  return _domain.getMin() +
         Eigen::Array3d((i + 0.5) * h[0], (j + 0.5) * h[1], (k + 0.5) * h[2]);
}

Eigen::Array3d CellCenteredGrid3::getCellPosition(int id) {
  auto ijk = idToijk(id);
  return getCellPosition(ijk[0], ijk[1], ijk[2]);
}

BoundingBox3 CellCenteredGrid3::getCellBoundingBox(int id) {
  auto ijk = idToijk(id);
  return getCellBoundingBox(ijk[0], ijk[1], ijk[2]);
}

BoundingBox3 CellCenteredGrid3::getCellBoundingBox(int i, int j, int k) {
  auto position = getCellPosition(i, j, k);
  BoundingBox3 box;
  Eigen::Array3d h = getH();
  box.setMin(position - h / 2);
  box.setMax(position + h / 2);
  return box;
}

std::vector<int> CellCenteredGrid3::idToijk(size_t id) {
  std::vector<int> indices(3);
  indices[2] = id / (_gridSize[0] * _gridSize[1]);
  indices[1] = (id % (_gridSize[0] * _gridSize[1])) / _gridSize[0];
  indices[0] = (id % (_gridSize[0] * _gridSize[1])) % _gridSize[0];
  return indices;
}

size_t CellCenteredGrid3::newCellScalarLabel(std::string label) {
  return newCellScalarLabel(label, 0);
}

size_t CellCenteredGrid3::newCellScalarLabel(std::string label, double value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _scalarData.emplace_back(std::vector<double>(_gridSize.prod(), value));
    _dataLabel[label] = _scalarData.size() - 1;
  }
  return _dataLabel[label];
}

size_t CellCenteredGrid3::newCellArrayLabel(std::string label) {
  return newCellArrayLabel(label, Eigen::Array3d(0, 0, 0));
}

size_t CellCenteredGrid3::newCellArrayLabel(std::string label,
                                            Eigen::Array3d value) {
  if (_dataLabel.find(label) == _dataLabel.end()) {
    _arrayData.emplace_back(
        std::vector<Eigen::Array3d>(_gridSize.prod(), value));
    _dataLabel[label] = _arrayData.size() - 1;
  }
  return _dataLabel[label];
}

std::vector<double> &CellCenteredGrid3::getCellScalarData(std::string label) {
  if (_dataLabel.find(label) != _dataLabel.end())
    return _scalarData[_dataLabel[label]];
  std::cerr << "\033[1;31m[ ERROR ]\033[0mLabel not found!\n";
  return _scalarData[0];
}

std::vector<double> &CellCenteredGrid3::getCellScalarData(size_t index) {
  return _scalarData[index];
}

std::vector<Eigen::Array3d> &
CellCenteredGrid3::getCellArrayData(std::string label) {
  return _arrayData[_dataLabel[label]];
}

std::vector<Eigen::Array3d> &CellCenteredGrid3::getCellArrayData(size_t index) {
  return _arrayData[index];
}

size_t CellCenteredGrid3::cellCount() { return _gridSize.prod(); }

double CellCenteredGrid3::interpolateCellScalarData(int dataId,
                                                    Eigen::Array3d position) {
  auto &data = getCellScalarData(dataId);
  auto h = getH();

  // Find which cell-id data belongs
  Eigen::Array3i cellId =
      (position - _domain.getMin()).cwiseQuotient(h).floor().cast<int>();

  // Assemble bilinear stencil interpolation for velocities
  auto cellPos = getCellPosition(cellId[0], cellId[1], cellId[2]);
  int index[3];
  for (size_t i = 0; i < 3; i++) {
    index[i] = std::min(_gridSize[i] - 1, cellId[i] + 1);
    if (position[i] < cellPos[i]) {
      index[i] = std::max(0, cellId[i] - 1);
      h[i] = -h[i];
    }
  }

  std::vector<Eigen::Array3d> points;
  std::vector<double> values(8);

  // cell Stencil for linear interpolation
  points.emplace_back(cellPos);
  points.emplace_back(cellPos + Eigen::Array3d(h[0], 0, 0));
  points.emplace_back(cellPos + Eigen::Array3d(0, h[1], 0));
  points.emplace_back(cellPos + Eigen::Array3d(h[0], h[1], 0));
  points.emplace_back(cellPos + Eigen::Array3d(0, 0, h[2]));
  points.emplace_back(cellPos + Eigen::Array3d(h[0], 0, h[2]));
  points.emplace_back(cellPos + Eigen::Array3d(0, h[1], h[2]));
  points.emplace_back(cellPos + h);

  values[0] = data[ijkToid(cellId[0], cellId[1], cellId[2])];
  values[1] = data[ijkToid(index[0], cellId[1], cellId[2])];
  values[2] = data[ijkToid(cellId[0], index[1], cellId[2])];
  values[3] = data[ijkToid(index[0], index[1], cellId[2])];
  values[4] = data[ijkToid(cellId[0], cellId[1], index[2])];
  values[5] = data[ijkToid(index[0], cellId[1], index[2])];
  values[6] = data[ijkToid(cellId[0], index[1], index[2])];
  values[7] = data[ijkToid(index[0], index[1], index[2])];

  return Ramuh::Interpolator::trilinear(position, points, values);
}

Eigen::Array3d
CellCenteredGrid3::interpolateCellArrayData(int dataId,
                                            Eigen::Array3d position) {
  auto &data = getCellArrayData(dataId);
  auto h = getH();

  // Find which cell-id data belongs
  Eigen::Array3i cellId =
      (position - _domain.getMin()).cwiseQuotient(h).floor().cast<int>();

  // Assemble bilinear stencil interpolation for velocities
  auto cellPos = getCellPosition(cellId[0], cellId[1], cellId[2]);
  int index[3];
  for (size_t i = 0; i < 3; i++) {
    index[i] = std::min(_gridSize[i] - 1, cellId[i] + 1);
    if (position[i] < cellPos[i]) {
      index[i] = std::max(0, cellId[i] - 1);
      h[i] = -h[i];
    }
  }

  std::vector<Eigen::Array3d> points;
  std::vector<double> values(4);

  // cell Stencil for linear interpolation
  points.emplace_back(cellPos);
  points.emplace_back(cellPos + Eigen::Array3d(h[0], 0, 0));
  points.emplace_back(cellPos + Eigen::Array3d(0, h[1], 0));
  points.emplace_back(cellPos + Eigen::Array3d(h[0], h[1], 0));
  points.emplace_back(cellPos + Eigen::Array3d(0, 0, h[2]));
  points.emplace_back(cellPos + Eigen::Array3d(h[0], 0, h[2]));
  points.emplace_back(cellPos + Eigen::Array3d(0, h[1], h[2]));
  points.emplace_back(cellPos + h);

  Eigen::Array3d interpData;
  for (size_t d = 0; d < 2; d++) {
    values[0] = data[ijkToid(cellId[0], cellId[1], cellId[2])][d];
    values[1] = data[ijkToid(index[0], cellId[1], cellId[2])][d];
    values[2] = data[ijkToid(cellId[0], index[1], cellId[2])][d];
    values[3] = data[ijkToid(index[0], index[1], cellId[2])][d];
    values[4] = data[ijkToid(cellId[0], cellId[1], index[2])][d];
    values[5] = data[ijkToid(index[0], cellId[1], index[2])][d];
    values[6] = data[ijkToid(cellId[0], index[1], index[2])][d];
    values[7] = data[ijkToid(index[0], index[1], index[2])][d];

    interpData[d] = Ramuh::Interpolator::trilinear(position, points, values);
  }

  return interpData;
}

} // namespace Ramuh