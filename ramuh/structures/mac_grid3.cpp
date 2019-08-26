#include <structures/mac_grid3.h>

namespace Ramuh {
MacGrid3::MacGrid3() : CellCenteredGrid3() {}

MacGrid3::MacGrid3(BoundingBox3 domain, Eigen::Array3i gridSize)
    : CellCenteredGrid3(domain, gridSize) {}

size_t MacGrid3::newFaceScalarLabel(std::string label) {
  return newFaceScalarLabel(label, 0);
}

size_t MacGrid3::newFaceScalarLabel(std::string label, double value) {
  if (_faceDataLabel.find(label) == _faceDataLabel.end()) {
    _uScalar.emplace_back(std::vector<double>(
        (_gridSize[0] + 1) * _gridSize[1] * _gridSize[2], value));
    _vScalar.emplace_back(std::vector<double>(
        _gridSize[0] * (_gridSize[1] + 1) * _gridSize[2], value));
    _wScalar.emplace_back(std::vector<double>(
        _gridSize[0] * _gridSize[1] * (_gridSize[2] + 1), value));
    _faceDataLabel[label] = _uScalar.size() - 1;
  }
  return _faceDataLabel[label];
}

std::vector<double> &MacGrid3::getFaceScalarData(size_t face,
                                                 std::string label) {
  if (face == 0)
    return _uScalar[_faceDataLabel[label]];
  if (face == 1)
    return _vScalar[_faceDataLabel[label]];
  return _wScalar[_faceDataLabel[label]];
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
  if (_faceDataLabel.find(label) == _faceDataLabel.end()) {
    _uArray.emplace_back(std::vector<Eigen::Array3d>(
        (_gridSize[0] + 1) * _gridSize[1] * _gridSize[2], value));
    _vArray.emplace_back(std::vector<Eigen::Array3d>(
        _gridSize[0] * (_gridSize[1] + 1) * _gridSize[2], value));
    _wArray.emplace_back(std::vector<Eigen::Array3d>(
        _gridSize[0] * _gridSize[1] * (_gridSize[2] + 1), value));
    _faceDataLabel[label] = _uArray.size() - 1;
  }
  return _faceDataLabel[label];
}

std::vector<Eigen::Array3d> &MacGrid3::getFaceArrayData(size_t face,
                                                        std::string label) {
  if (face == 0)
    return _uArray[_faceDataLabel[label]];
  if (face == 1)
    return _vArray[_faceDataLabel[label]];
  return _wArray[_faceDataLabel[label]];
}

std::vector<Eigen::Array3d> &MacGrid3::getFaceArrayData(size_t face,
                                                        size_t index) {
  if (face == 0)
    return _uArray[index];
  if (face == 1)
    return _vArray[index];
  return _wArray[index];
}

Eigen::Array3d MacGrid3::facePosition(int face, int id) {
  auto ijk = faceIdToijk(face, id);
  auto h = getH();
  int i = ijk[0], j = ijk[1], k = ijk[2];
  if (face == 2)
    return _domain.min() + Eigen::Array3d(i * h[0], j * h[1], (k + 0.5) * h[2]);
  if (face == 1)
    return _domain.min() + Eigen::Array3d(i * h[0], (j + 0.5) * h[1], k * h[2]);
  if (face == 0)
    return _domain.min() + Eigen::Array3d((i + 0.5) * h[0], j * h[1], k * h[2]);
  return Eigen::Array3d(0);
}

Eigen::Array3d MacGrid3::facePosition(int face, int i, int j, int k) {
  return facePosition(face, faceijkToid(face, i, j, k));
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

} // namespace Ramuh
