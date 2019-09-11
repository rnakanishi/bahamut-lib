
#include <blas/weno.h>
#include <structures/mac_grid2.h>
#include <cmath>
#include <iostream>

class DiffGrid : public Ramuh::MacGrid2 {
public:
  DiffGrid()
      : Ramuh::MacGrid2(Ramuh::BoundingBox2(Eigen::Array2d(0, 0),
                                            Eigen::Array2d(M_PI, M_PI)),
                        Eigen::Array2i(10, 10)) {}

  void initialize() {
    _functionId = newScalarLabel("function");
    _analyticId = newScalarLabel("analytic");
    _wenoId = newScalarLabel("weno");

    _velocityId = newFaceArrayLabel("velocity");

    auto &function = getScalarLabel(_functionId);
    auto &analytic = getScalarLabel(_analyticId);
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getPosition(i);
      function[i] = sin(p[0]) * cos(p[1]);
      analytic[i] = (cos(p[0]) * cos(p[1]));
    }

    auto &velocity = getFaceArrayLabel(0, _velocityId);
    for (size_t i = 0; i < uFaceCount(); i++) {
      velocity[i] = 1.0;
    }
  }

  void solveODE() {
    auto h = getH();
    auto &phi = getScalarLabel(_functionId);
    auto &analytic = getScalarLabel(_analyticId);
    auto &weno = getScalarLabel(_wenoId);

    // Weno computation
    std::vector<double> values(6);
    for (int id = 0; id < cellCount(); id++) {
      auto p = getPosition(id);
      auto ij = idToij(id);
      int i = ij.first, j = ij.second;

      for (int ival = 0, ii = 3; ival < 7; ival++, ii--) {
        int index = std::min(_gridSize[0] - 1, std::max(1, i + ii));
        int index1 = std::min(_gridSize[0] - 2, std::max(0, i + ii - 1));
        values[ival] = (phi[ijToid(index, j)] - phi[ijToid(index1, j)]) / h[0];
      }
      if (i == 0 || i == 1)
        values[5] = values[4] = 1e4;
      if (i == _gridSize[0] - 1 || i == _gridSize[0] - 2)
        values[0] = values[1] = 1e4;
      double dPhi = Ramuh::Weno::evaluate(values, h[0], false);
      weno[id] = dPhi;
    }
    for (size_t i = 0; i < cellCount(); i++) {
      std::cout << analytic[i] << " " << weno[i] << "; ";
    }
    std::cout << std::endl;
  }

private:
  size_t _functionId;
  size_t _analyticId;
  size_t _wenoId;
  size_t _velocityId;
};

int main(int argc, char const *argv[]) {
  DiffGrid grid;
  grid.initialize();
  grid.solveODE();

  return 0;
}
