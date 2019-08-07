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

std::vector<double> &MacGrid3::getFaceScalarVector(size_t face,
                                                   std::string label) {
  if (face == 0)
    return _uScalar[_faceDataLabel[label]];
  if (face == 1)
    return _vScalar[_faceDataLabel[label]];
  return _wScalar[_faceDataLabel[label]];
}

std::vector<double> &MacGrid3::getFaceScalarVector(size_t face, size_t index) {
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

std::vector<Eigen::Array3d> &MacGrid3::getFaceArrayVector(size_t face,
                                                          std::string label) {
  if (face == 0)
    return _uArray[_faceDataLabel[label]];
  if (face == 1)
    return _vArray[_faceDataLabel[label]];
  return _wArray[_faceDataLabel[label]];
}

std::vector<Eigen::Array3d> &MacGrid3::getFaceArrayVector(size_t face,
                                                          size_t index) {
  if (face == 0)
    return _uArray[index];
  if (face == 1)
    return _vArray[index];
  return _wArray[index];
}

size_t MacGrid3::cellCount() { return _gridSize.prod(); }

} // namespace Ramuh
