#include <structures/mac_grid2.h>

namespace Ramuh {
MacGrid2::MacGrid2() : CellCenteredGrid2() {}

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

std::vector<double> &MacGrid2::getFaceScalarLabel(size_t face,
                                                  std::string label) {
  if (face == 0)
    return _uScalar[_faceDataLabel[label]];
  return _vScalar[_faceDataLabel[label]];
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

std::vector<Eigen::Array2d> &MacGrid2::getFaceArrayLabel(size_t face,
                                                         std::string label) {
  if (face == 0)
    return _uArray[_faceDataLabel[label]];
  return _vArray[_faceDataLabel[label]];
}

} // namespace Ramuh
