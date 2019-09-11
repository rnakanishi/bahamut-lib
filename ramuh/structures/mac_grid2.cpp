#include <structures/mac_grid2.h>

namespace Ramuh {
MacGrid2::MacGrid2() : CellCenteredGrid2() {}

MacGrid2::MacGrid2(BoundingBox2 domain, Eigen::Array2i gridSize)
    : CellCenteredGrid2(domain, gridSize) {}

size_t MacGrid2::newFaceScalarLabel(std::string label) {
  return newFaceScalarLabel(label, 0);
}

size_t MacGrid2::newFaceScalarLabel(std::string label, double value) {
  if (_faceDataLabel.find(label) == _faceDataLabel.end()) {
    _uScalar.emplace_back(
        std::vector<double>((_gridSize[0] + 1) * _gridSize[1], value));
    _vScalar.emplace_back(
        std::vector<double>(_gridSize[0] * (_gridSize[1] + 1), value));
    _faceDataLabel[label] = _uScalar.size() - 1;
  }
  return _faceDataLabel[label];
}

std::vector<double> &MacGrid2::getFaceScalarData(size_t face,
                                                 std::string label) {
  if (face == 0)
    return _uScalar[_faceDataLabel[label]];
  return _vScalar[_faceDataLabel[label]];
}

std::vector<double> &MacGrid2::getFaceScalarData(size_t face, int id) {
  if (face == 0)
    return _uScalar[id];
  return _vScalar[id];
}

size_t MacGrid2::newFaceArrayLabel(std::string label) {
  return newFaceArrayLabel(label, Eigen::Array2d(0., 0.));
}

size_t MacGrid2::newFaceArrayLabel(std::string label, Eigen::Array2d value) {
  if (_faceDataLabel.find(label) == _faceDataLabel.end()) {
    _uArray.emplace_back(
        std::vector<Eigen::Array2d>((_gridSize[0] + 1) * _gridSize[1], value));
    _vArray.emplace_back(
        std::vector<Eigen::Array2d>(_gridSize[0] * (_gridSize[1] + 1), value));
    _faceDataLabel[label] = _uArray.size() - 1;
  }
  return _faceDataLabel[label];
}

std::vector<Eigen::Array2d> &MacGrid2::getFaceArrayData(size_t face,
                                                        std::string label) {
  if (face == 0)
    return _uArray[_faceDataLabel[label]];
  return _vArray[_faceDataLabel[label]];
}

std::vector<Eigen::Array2d> &MacGrid2::getFaceArrayData(size_t face, int id) {
  if (face == 0)
    return _uArray[id];
  return _vArray[id];
}

int MacGrid2::faceCount(int face) {
  if (face == 0)
    return (_gridSize[0] + 1) * _gridSize[1];
  else
    return _gridSize[0] * (_gridSize[1] + 1);
}

Eigen::Array2d MacGrid2::facePosition(size_t face, int i, int j) {
  return facePosition(face, faceijToid(face, i, j));
}

Eigen::Array2d MacGrid2::facePosition(size_t face, int faceId) {
  auto ij = faceIdToij(face, faceId);
  auto h = getH();
  int i = ij[0], j = ij[1];
  if (face == 0)
    return _domain.min() + Eigen::Array2d((i + 0.5) * h[0], j * h[1]);
  return _domain.min() + Eigen::Array2d(i * h[0], (j + 0.5) * h[1]);
}

int MacGrid2::faceijToid(int face, int i, int j) {
  if (face == 0)
    return j * (_gridSize[0] + 1) + i;
  else
    return j * _gridSize[0] + i;
}

std::vector<int> MacGrid2::faceIdToij(int face, int id) {
  std::vector<int> index(2);
  if (face == 0) {
    index[1] = id / (_gridSize[0] + 1);
    index[0] = id % (_gridSize[0] + 1);
  } else {
    index[1] = id / (_gridSize[0]);
    index[0] = id % (_gridSize[0]);
  }
  return index;
}

} // namespace Ramuh
