#include <structures/mac_grid2.h>

namespace Ramuh {
MacGrid2::MacGrid2() : CellCenteredGrid2() {}

size_t MacGrid2::newFaceLabel(std::string label) {
  return newFaceLabel(label, 0);
}

size_t MacGrid2::newFaceLabel(std::string label, double value) {
  if (_faceDataLabel.find(label) == _faceDataLabel.end()) {
    _udata.emplace_back(
        std::vector<double>((_gridSize[0] + 1) * _gridSize[1], value));
    _vdata.emplace_back(
        std::vector<double>(_gridSize[0] * (_gridSize[1] + 1), value));
    _faceDataLabel[label] = _udata.size() - 1;
  }
  return _faceDataLabel[label];
}

std::vector<double> &MacGrid2::getFaceLabel(size_t face, std::string label) {
  if (face == 0)
    return _udata[_faceDataLabel[label]];
  return _vdata[_faceDataLabel[label]];
}

} // namespace Ramuh
