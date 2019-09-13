#include <blas/weno.h>
#include <structures/mac_grid2.h>
#include <structures/mac_grid1.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

class DiffGrid : public Ramuh::MacGrid1 {
public:
  DiffGrid() : Ramuh::MacGrid1(Ramuh::BoundingBox1(0, 2 * M_PI), 50) {}

  void initialize() {
    _dt = 1 / 60.;
    _ellapsed = 0.;

    _functionId = newLabel("function");
    _analyticId = newLabel("analyticSolution");
    _discreteId = newLabel("wenoSolution");

    auto &values = getScalarData(_functionId);
    auto &analytic = getScalarData(_analyticId);
    auto h = getH();
    for (int i = 0; i < _gridSize; i++) {
      double position = getPosition(i);
      values[i] = std::cos(position);
      analytic[i] = cos(position);
    }
  }

  void solveTimestep() {
    auto &phi = getScalarData(_functionId);
    auto &analytic = getScalarData(_analyticId);
    auto &weno = getScalarData(_discreteId);
    auto h = getH();
    std::vector<double> values(6);

    for (int i = 0; i < _gridSize; i++) {
      double position = getPosition(i);
      double x = position;

      analytic[i] = cos(position - _ellapsed);
      // compute weno
      // Careful with velocity orientation and ii sign
      for (int ival = 0, ii = -3; ival < 7; ival++, ii++) {
        int index = std::min(_gridSize - 1, std::max(1, i + ii + 1));
        int index1 = std::min(_gridSize - 2, std::max(0, i + ii));
        values[ival] = (phi[index] - phi[index1]) / h;
      }
      //   boundary conditions
      if (i == 0 || i == 1)
        values[5] = values[4] = 1e6;
      if (i == _gridSize - 1 || i == _gridSize - 2)
        values[0] = values[1] = 1e6;

      double dPhi = Ramuh::Weno::evaluate(values, h, false);
      if (i == _gridSize - 1)
        dPhi = Ramuh::Weno::evaluate(
            std::vector<double>(values.begin() + 1, values.end()), h, false);
      weno[i] = dPhi;
    }

    // Time integration
    for (int i = 0; i < _gridSize; i++) {
      phi[i] = phi[i] - weno[i] * _dt;
    }
    _ellapsed += _dt;
  }

  void print() {
    auto &analytic = getScalarData(_analyticId);
    auto &phi = getScalarData(_functionId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/weno/1d/" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);

    for (int i = 0; i < _gridSize; i++) {
      file << analytic[i] << " " << phi[i] << ";\n";
    }
    file.close();
  }

private:
  int _functionId;
  int _analyticId, _discreteId;
  double _dt, _ellapsed;
};

int main(int argc, char const *argv[]) {

  DiffGrid grid;
  grid.initialize();
  for (size_t i = 0; i < 90; i++) {
    grid.solveTimestep();
    grid.print();
  }

  return 0;
}
