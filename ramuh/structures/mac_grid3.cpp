#include <structures/mac_grid3.h>
#include <blas/interpolator.h>

namespace Ramuh {
MacGrid3::MacGrid3() : CellCenteredGrid3() {}

MacGrid3::MacGrid3(BoundingBox3 domain, Eigen::Array3i gridSize)
    : CellCenteredGrid3(domain, gridSize) {}

size_t MacGrid3::newFaceScalarLabel(std::string label) {
  return newFaceScalarLabel(label, 0);
}

size_t MacGrid3::newFaceScalarLabel(std::string label, double value) {
  if (_faceScalarLabel.find(label) == _faceScalarLabel.end()) {
    _uScalar.emplace_back(std::vector<double>(
        (_gridSize[0] + 1) * _gridSize[1] * _gridSize[2], value));
    _vScalar.emplace_back(std::vector<double>(
        _gridSize[0] * (_gridSize[1] + 1) * _gridSize[2], value));
    _wScalar.emplace_back(std::vector<double>(
        _gridSize[0] * _gridSize[1] * (_gridSize[2] + 1), value));
    _faceScalarLabel[label] = _uScalar.size() - 1;
  }
  return _faceScalarLabel[label];
}

std::vector<double> &MacGrid3::getFaceScalarData(size_t face,
                                                 std::string label) {
  if (face == 0)
    return _uScalar[_faceScalarLabel[label]];
  if (face == 1)
    return _vScalar[_faceScalarLabel[label]];
  return _wScalar[_faceScalarLabel[label]];
}

std::vector<double> &MacGrid3::getFaceScalarData(size_t face, size_t index) {
  if (face == 0)
    return _uScalar[index];
  if (face == 1)
    return _vScalar[index];
  return _wScalar[index];
}

size_t MacGrid3::newFaceArrayLabel(std::string label) {
  return newFaceArrayLabel(label, Eigen::Array3d(0., 0., 0.));
}

size_t MacGrid3::newFaceArrayLabel(std::string label, Eigen::Array3d value) {
  if (_faceArrayLabel.find(label) == _faceArrayLabel.end()) {
    _uArray.emplace_back(std::vector<Eigen::Array3d>(
        (_gridSize[0] + 1) * _gridSize[1] * _gridSize[2], value));
    _vArray.emplace_back(std::vector<Eigen::Array3d>(
        _gridSize[0] * (_gridSize[1] + 1) * _gridSize[2], value));
    _wArray.emplace_back(std::vector<Eigen::Array3d>(
        _gridSize[0] * _gridSize[1] * (_gridSize[2] + 1), value));
    _faceArrayLabel[label] = _uArray.size() - 1;
  }
  return _faceArrayLabel[label];
}

int MacGrid3::getFaceScalarLabelId(std::string label) {
  if (_faceScalarLabel.find(label) == _faceScalarLabel.end())
    return -1;
  return _faceScalarLabel[label];
}

int MacGrid3::getFaceArrayLabelId(std::string label) {
  if (_faceArrayLabel.find(label) == _faceArrayLabel.end())
    return -1;
  return _faceArrayLabel[label];
}

std::vector<Eigen::Array3d> &MacGrid3::getFaceArrayData(size_t face,
                                                        std::string label) {
  if (face == 0)
    return _uArray[_faceArrayLabel[label]];
  if (face == 1)
    return _vArray[_faceArrayLabel[label]];
  return _wArray[_faceArrayLabel[label]];
}

std::vector<Eigen::Array3d> &MacGrid3::getFaceArrayData(size_t face,
                                                        size_t index) {
  if (face == 0)
    return _uArray[index];
  if (face == 1)
    return _vArray[index];
  return _wArray[index];
}

Eigen::Array3d MacGrid3::getFacePosition(int face, int id) {
  auto ijk = faceIdToijk(face, id);
  auto h = getH();
  int i = ijk[0], j = ijk[1], k = ijk[2];
  if (face == 2)
    return _domain.getMin() +
           Eigen::Array3d((i + .5) * h[0], (j + .5) * h[1], k * h[2]);
  if (face == 1)
    return _domain.getMin() +
           Eigen::Array3d((i + .5) * h[0], j * h[1], (k + .5) * h[2]);
  if (face == 0)
    return _domain.getMin() +
           Eigen::Array3d(i * h[0], (j + .5) * h[1], (k + .5) * h[2]);
  return Eigen::Array3d(0);
}

Eigen::Array3d MacGrid3::getFacePosition(int face, int i, int j, int k) {
  return getFacePosition(face, faceijkToid(face, i, j, k));
}

int MacGrid3::faceCount(int face) {
  if (face == 2)
    return _gridSize[0] * _gridSize[1] * (_gridSize[2] + 1);
  else if (face == 1)
    return _gridSize[0] * (_gridSize[1] + 1) * _gridSize[2];
  else
    return (_gridSize[0] + 1) * _gridSize[1] * _gridSize[2];
}

int MacGrid3::faceijkToid(int face, int i, int j, int k) {
  if (face == 2)
    return k * _gridSize[0] * _gridSize[1] + j * _gridSize[0] + i;
  else if (face == 1)
    return k * _gridSize[0] * (_gridSize[1] + 1) + j * _gridSize[0] + i;
  else
    return k * (_gridSize[0] + 1) * _gridSize[1] + j * (_gridSize[0] + 1) + i;
}

std::vector<int> MacGrid3::faceIdToijk(int face, int id) {
  std::vector<int> indices(3);
  if (face == 2) {
    indices[2] = id / (_gridSize[0] * _gridSize[1]);
    indices[1] = (id % (_gridSize[0] * _gridSize[1])) / _gridSize[0];
    indices[0] = (id % (_gridSize[0] * _gridSize[1])) % _gridSize[0];
  } else if (face == 1) {
    indices[2] = id / (_gridSize[0] * (_gridSize[1] + 1));
    indices[1] = (id % (_gridSize[0] * (_gridSize[1] + 1))) / _gridSize[0];
    indices[0] = (id % (_gridSize[0] * (_gridSize[1] + 1))) % _gridSize[0];
  } else {
    indices[2] = id / ((_gridSize[0] + 1) * _gridSize[1]);
    indices[1] =
        (id % ((_gridSize[0] + 1) * _gridSize[1])) / (_gridSize[0] + 1);
    indices[0] =
        (id % ((_gridSize[0] + 1) * _gridSize[1])) % (_gridSize[0] + 1);
  }
  return indices;
}

Eigen::Array3i MacGrid3::getFaceGridSize(int face) {
  Eigen::Array3i inc(0);
  inc[face] += 1;
  return _gridSize + inc;
}

double MacGrid3::interpolateFaceScalarData(int face, int dataId,
                                           Eigen::Array3d position) {
  auto &data = getFaceScalarData(face, dataId);
  auto h = getH();

  Eigen::Array3i cellId =
      (position - _domain.getMin()).cwiseQuotient(h).floor().cast<int>();

  // Assemble bilinear stencil interpolation for velocities
  auto facePos = getFacePosition(face, cellId[0], cellId[1], cellId[2]);
  auto faceGridSize = getFaceGridSize(face);
  int index[3];
  for (size_t i = 0; i < 3; i++) {
    index[i] = std::min(faceGridSize[i] - 1, cellId[i] + 1);
    if (position[i] < facePos[i]) {
      index[i] = cellId[i];
      cellId[i] = std::max(0, cellId[i] - 1);
      facePos[i] -= h[i];
    }
  }

  std::vector<Eigen::Array3d> points(8);
  std::vector<double> values(8);

  points[0] = (facePos + Eigen::Array3d(0, 0, 0));
  points[1] = (facePos + Eigen::Array3d(h[0], 0, 0));
  points[2] = (facePos + Eigen::Array3d(0, h[1], 0));
  points[3] = (facePos + Eigen::Array3d(h[0], h[1], 0));
  points[4] = (facePos + Eigen::Array3d(0, 0, h[2]));
  points[5] = (facePos + Eigen::Array3d(h[0], 0, h[2]));
  points[6] = (facePos + Eigen::Array3d(0, h[1], h[2]));
  points[7] = (facePos + Eigen::Array3d(h[0], h[1], h[2]));

  values[0] = data[faceijkToid(face, cellId[0], cellId[1], cellId[2])];
  values[1] = data[faceijkToid(face, index[0], cellId[1], cellId[2])];
  values[2] = data[faceijkToid(face, cellId[0], index[1], cellId[2])];
  values[3] = data[faceijkToid(face, index[0], index[1], cellId[2])];
  values[4] = data[faceijkToid(face, cellId[0], cellId[1], index[2])];
  values[5] = data[faceijkToid(face, index[0], cellId[1], index[2])];
  values[6] = data[faceijkToid(face, cellId[0], index[1], index[2])];
  values[7] = data[faceijkToid(face, index[0], index[1], index[2])];

  return Ramuh::Interpolator::trilinear(position, points, values);
}

Eigen::Array3d MacGrid3::interpolateFaceArrayData(int face, int dataId,
                                                  Eigen::Array3d position) {
  auto &data = getFaceArrayData(face, dataId);
  auto h = getH();

  Eigen::Array3i cellId =
      (position - _domain.getMin()).cwiseQuotient(h).floor().cast<int>();

  // Assemble bilinear stencil interpolation for velocities
  auto facePos = getFacePosition(face, cellId[0], cellId[1], cellId[2]);
  auto faceGridSize = getFaceGridSize(face);
  int index[3];
  for (size_t i = 0; i < 3; i++) {
    index[i] = std::min(faceGridSize[i] - 1, cellId[i] + 1);
    if (position[i] < facePos[i]) {
      index[i] = std::max(0, cellId[i] - 1);
      h[i] = -h[i];
    }
  }

  std::vector<Eigen::Array3d> points;
  std::vector<double> values(8);

  points.emplace_back(facePos);
  points.emplace_back(facePos + Eigen::Array3d(h[0], 0, 0));
  points.emplace_back(facePos + Eigen::Array3d(0, h[1], 0));
  points.emplace_back(facePos + Eigen::Array3d(h[0], h[1], 0));
  points.emplace_back(facePos + Eigen::Array3d(0, 0, h[2]));
  points.emplace_back(facePos + Eigen::Array3d(h[0], 0, h[2]));
  points.emplace_back(facePos + Eigen::Array3d(0, h[1], h[2]));
  points.emplace_back(facePos + h);

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
