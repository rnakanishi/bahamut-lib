#include <structures/mac_grid1.h>

namespace Ramuh {
MacGrid1::MacGrid1() : CellCenteredGrid1() {}

MacGrid1::MacGrid1(BoundingBox1 domain, int gridSize)
    : CellCenteredGrid1(domain, gridSize) {}

size_t MacGrid1::newFaceLabel(std::string label) {
  return newFaceLabel(label, 0);
}

size_t MacGrid1::newFaceLabel(std::string label, double value) {
  if (_faceLabel.find(label) == _faceLabel.end()) {
    _uData.emplace_back(std::vector<double>(_gridSize + 1, value));
    _faceLabel[label] = _uData.size() - 1;
  }
  return _faceLabel[label];
}

std::vector<double> &MacGrid1::getFaceLabelData(std::string label) {
  return _uData[_faceLabel[label]];
}

std::vector<double> &MacGrid1::getFaceLabelData(int id) { return _uData[id]; }

} // namespace Ramuh
