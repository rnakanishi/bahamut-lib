#include <blas/weno.h>
#include <structures/mac_grid2.h>
#include <structures/mac_grid1.h>
#include <cmath>
#include <iostream>
#include <utility>
#include <cmath>

class DiffGrid : public Ramuh::MacGrid1 {
public:
  DiffGrid() : Ramuh::MacGrid1(Ramuh::BoundingBox1(0, M_PI), 20) {}

  void initialize() {
    _functionId = newLabel("function");
    _analyticId = newLabel("analyticSolution");
    _discreteId = newLabel("wenoSolution");
  }

  void setInitialValues() {
    auto &values = getLabelData(_functionId);
    auto h = getH();
    for (int i = 0; i < _gridSize; i++) {
      double position = getPosition(i);
      values[i] = std::sin(position);
    }
  }

  void compareSolution() {
    auto &phi = getLabelData(_functionId);
    auto &analytic = getLabelData(_analyticId);
    auto &weno = getLabelData(_discreteId);
    auto h = getH();
    std::vector<double> values(6);

    for (int i = 0; i < _gridSize; i++) {
      auto h = getH();
      double position = getPosition(i);
      double x = position;
      double solution = -cos(position);
      analytic[i] = solution;

      // compute weno
      for (int ival = 0, ii = 3; ival < 7; ival++, ii--) {
        int index = std::min(_gridSize - 1, std::max(0, i + ii));
        int index1 = std::min(_gridSize - 1, std::max(0, i + ii - 1));
        values[ival] = (phi[index] + phi[index1]) / 2;
      }
      double dPhi = Ramuh::Weno::evaluate(values, h, false);
      weno[i] = dPhi;
    }
    int newWenoId = newLabel("newWeno");
    auto &newWeno = getLabelData(newWenoId);

    double alpha = 0;
    for (size_t i = 0; i < _gridSize - 1; i++) {
      if (i > 0)
        alpha = std::max(alpha, abs(weno[i] - weno[i - 1]) / h);
      else
        alpha = std::max(alpha, abs(weno[i + 1] - weno[i]) / h);
    }

    for (size_t i = 0; i < _gridSize - 1; i++) {
      auto dPhi = weno[i];
      auto dPhi1 = weno[i + 1];
      // Monotone flux - Lax Friedrichs
      double flux = 0.5 * (dPhi + dPhi1);
      newWeno[i] = flux;
    }
    for (size_t i = 0; i < _gridSize - 1; i++) {
      if (i > 0)
        weno[i] = -(newWeno[i] - newWeno[i - 1]) / h;
      else
        weno[i] = -(newWeno[i + 1] - newWeno[i]) / h;
    }

    for (int i = 0; i < _gridSize; i++) {
      std::cout << phi[i] << " " << analytic[i] << " " << weno[i] << std::endl;
    }
  }

private:
  int _functionId;
  int _analyticId, _discreteId;
};

int main(int argc, char const *argv[]) {

  DiffGrid grid;
  grid.initialize();
  grid.setInitialValues();
  grid.compareSolution();

  return 0;
}
