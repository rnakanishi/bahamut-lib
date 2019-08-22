
#include <blas/weno.h>
#include <structures/mac_grid2.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

class DiffGrid : public Ramuh::MacGrid2 {
public:
  DiffGrid()
      : Ramuh::MacGrid2(Ramuh::BoundingBox2(Eigen::Array2d(0, 0),
                                            Eigen::Array2d(2 * M_PI, 2 * M_PI)),
                        Eigen::Array2i(25, 25)) {
    _dt = 1 / 60.;
    _ellapsed = 0.;
  }

  void initialize() {
    _functionId = newScalarLabel("function");
    _analyticId = newScalarLabel("analytic");
    _wenoId = newScalarLabel("weno");

    _velocityId = newFaceArrayLabel("velocity");

    auto &function = getScalarLabel(_functionId);
    auto &analytic = getScalarLabel(_analyticId);
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getPosition(i);
      function[i] = analytic[i] = sin(p[1]);
    }

    auto &velocity = getFaceArrayLabel(1, _velocityId);
    for (size_t i = 0; i < uFaceCount(); i++) {
      velocity[i] = 1.0;
    }
  }

  void solveTimestep() {
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

      analytic[id] = sin(p[1] - _ellapsed);

      for (int ival = 0, ii = -3; ival < 7; ival++, ii++) {
        int index = std::min(_gridSize[1] - 1, std::max(1, j + ii + 1));
        int index1 = std::min(_gridSize[1] - 2, std::max(0, j + ii));
        values[ival] = (phi[ijToid(i, index)] - phi[ijToid(i, index1)]) / h[1];
      }
      if (j == 0 || j == 1)
        values[5] = values[4] = 1e4;
      if (j == _gridSize[1] - 1 || j == _gridSize[1] - 2)
        values[0] = values[1] = 1e4;
      double dPhi = Ramuh::Weno::evaluate(values, h[1], false);
      weno[id] = dPhi;
    }

    for (int id = 0; id < cellCount(); id++) {
      phi[id] = phi[id] - weno[id] * _dt;
    }
    _ellapsed += _dt;
  }

  void print() {
    auto &analytic = getScalarLabel(_analyticId);
    auto &phi = getScalarLabel(_functionId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/weno/2d/" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);

    for (size_t i = 0; i < cellCount(); i++) {
      file << analytic[i] << " " << phi[i] << ";\n ";
    }
    file.close();
  }

private:
  size_t _functionId;
  size_t _analyticId;
  size_t _wenoId;
  size_t _velocityId;

  double _dt, _ellapsed;
};

int main(int argc, char const *argv[]) {
  DiffGrid grid;
  grid.initialize();
  for (size_t i = 0; i < 500; i++) {
    grid.solveTimestep();
    grid.print();
  }

  return 0;
}
