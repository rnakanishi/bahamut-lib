#include <structures/mac_grid3.h>

namespace Ramuh {
MacGrid3::MacGrid3() : CellCenteredGrid3() {}

size_t MacGrid3::newFaceLabel(std::string label) {
  return newFaceLabel(label, 0);
}

size_t MacGrid3::newFaceLabel(std::string label, double value) {
  if (_faceDataLabel.find(label) == _faceDataLabel.end()) {
    _udata.emplace_back(std::vector<double>(
        (_gridSize[0] + 1) * _gridSize[1] * _gridSize[2], value));
    _vdata.emplace_back(std::vector<double>(
        _gridSize[0] * (_gridSize[1] + 1) * _gridSize[2], value));
    _wdata.emplace_back(std::vector<double>(
        _gridSize[0] * _gridSize[1] * (_gridSize[2] + 1), value));
    _faceDataLabel[label] = _udata.size() - 1;
  }
  return _faceDataLabel[label];
}

std::vector<double> &MacGrid3::getFaceLabel(size_t face, std::string label) {
  if (face == 0)
    return _udata[_faceDataLabel[label]];
  if (face == 1)
    return _vdata[_faceDataLabel[label]];
  return _wdata[_faceDataLabel[label]];
}

} // namespace Ramuh
