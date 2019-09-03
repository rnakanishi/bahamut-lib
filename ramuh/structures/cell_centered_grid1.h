#ifndef __RAMUH_CELL_CENTERED_GRID_1_H__
#define __RAMUH_CELL_CENTERED_GRID_1_H__

#include <vector>
#include <map>
#include <string>
#include <geometry/bounding_box.h>

namespace Ramuh {

class CellCenteredGrid1 {
public:
  CellCenteredGrid1();
  CellCenteredGrid1(BoundingBox1 domain, int gridSize);

  void setGridSize(int size);

  double getH();

  double getPosition(int i);

  int cellCount();

  size_t newLabel(std::string label);
  size_t newLabel(std::string label, double value);

  std::vector<double> &getScalarData(std::string label);
  std::vector<double> &getScalarData(size_t id);

protected:
  int _gridSize;
  BoundingBox1 _domain;

  std::vector<std::vector<double>> _scalarData;
  std::map<std::string, size_t> _dataLabel;
};

} // namespace Ramuh

#endif