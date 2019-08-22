
#include <blas/weno.h>
#include <structures/mac_grid2.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

class DiffGrid : public Ramuh::MacGrid2 {
public:
  DiffGrid()
      : Ramuh::MacGrid2(Ramuh::BoundingBox2(Eigen::Array2d(-M_PI, -M_PI),
                                            Eigen::Array2d(M_PI, M_PI)),
                        Eigen::Array2i(50, 50)) {
    _dt = 1 / 100.;
    _ellapsedTime = 0.;
  }

  void initialize() {
    _functionId = newScalarLabel("function");
    _analyticId = newScalarLabel("analytic");
    _wenoId = newScalarLabel("weno");

    _velocityId = newFaceScalarLabel("velocity");

    auto &function = getScalarData(_functionId);
    auto &analytic = getScalarData(_analyticId);
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getPosition(i).abs();

      //   if (p.maxCoeff() <= 1)
      function[i] = analytic[i] = p.maxCoeff() - 1;
      //   else
      //     function[i] = analytic[i] = 0;
    }

    auto &u = getFaceScalarData(0, _velocityId);
    for (size_t i = 0; i < faceCount(0); i++) {
      auto ij = faceIdToij(0, i);
      auto p = facePosition(0, i);
      if (ij.first > 0 && ij.first < _gridSize[0])
        u[i] = -p[1];
      //   else
      //     u[i] = 0;
    }
    auto &v = getFaceScalarData(1, _velocityId);
    for (size_t i = 0; i < faceCount(1); i++) {
      auto ij = faceIdToij(1, i);
      auto p = facePosition(1, i);
      if (ij.second > 0 && ij.second < _gridSize[1])
        v[i] = p[0];
      //   else
      //     v[i] = 0;
    }
  }

  void solveTimestep() {
    auto h = getH();
    auto &phi = getScalarData(_functionId);
    auto &analytic = getScalarData(_analyticId);
    auto &weno = getScalarData(_wenoId);

    // Weno computation
    std::vector<double> values(6);
    for (int id = 0; id < cellCount(); id++) {
      auto p = getPosition(id);
      auto ij = idToij(id);
      int i = ij.first, j = ij.second;
      weno[id] = 0;

      Eigen::Matrix2d rotation;
      double angle = _ellapsedTime;
      rotation.row(0) << cos(angle), sin(angle);
      rotation.row(1) << -sin(angle), cos(angle);
      p = rotation * p.matrix();
      //   analytic[id] = sin(p[0]) * cos(p[1]);
      //   analytic[id] = sin(p[0] + _ellapsedTime) * cos(p[1] + _ellapsedTime);

      for (size_t coord = 0; coord < 2; coord++) {
        auto &u = getFaceScalarData(coord, _velocityId);
        auto faceij = ij;
        int facei = faceij.first, facej = faceij.second;

        double velocity;
        if (coord == 0)
          velocity = ((u[faceijToid(coord, facei, facej)] +
                       u[faceijToid(coord, facei + 1, facej)]) /
                      2);
        else
          velocity = ((u[faceijToid(coord, facei, facej)] +
                       u[faceijToid(coord, facei, facej + 1)]) /
                      2);
        if (velocity <= 0) { // Downwind
          for (int ival = 0, inc = 3; ival < 7; ival++, inc--) {
            int index, index1, faceIndex;
            if (coord == 0) {
              index = std::min(_gridSize[coord] - 1, std::max(1, i + inc));
              index1 = std::min(_gridSize[coord] - 2, std::max(0, i + inc - 1));
              faceIndex = std::min(_gridSize[coord], std::max(0, i + inc));
              values[ival] = u[faceijToid(coord, faceIndex, facej)] *
                             (phi[ijToid(index, j)] - phi[ijToid(index1, j)]) /
                             h[coord];
              if (i <= 1)
                values[5] = 1e4;
              if (i <= 0)
                values[4] = 1e4;
              if (i >= _gridSize[coord] - 2)
                values[0] = 1e4;
              if (i >= _gridSize[coord] - 1)
                values[1] = 1e4;
            } else {
              index = std::min(_gridSize[coord] - 1, std::max(1, j + inc));
              index1 = std::min(_gridSize[coord] - 2, std::max(0, j + inc - 1));
              faceIndex = std::min(_gridSize[coord], std::max(0, j + inc));
              values[ival] = u[faceijToid(coord, facei, faceIndex)] *
                             (phi[ijToid(i, index)] - phi[ijToid(i, index1)]) /
                             h[coord];
              if (j <= 1)
                values[5] = 1e4;
              if (j <= 0)
                values[4] = 1e4;
              if (j >= _gridSize[coord] - 2)
                values[0] = 1e4;
              if (j >= _gridSize[coord] - 1)
                values[1] = 1e4;
            }
          }
        } else { // Upwind
          for (int ival = 0, inc = -3; ival < 7; ival++, inc++) {
            int index, index1, faceIndex;
            if (coord == 0) {
              index = std::min(_gridSize[coord] - 1, std::max(1, i + inc + 1));
              index1 = std::min(_gridSize[coord] - 2, std::max(0, i + inc));
              faceIndex = std::min(_gridSize[coord], std::max(0, i + inc + 1));
              values[ival] = u[faceijToid(coord, faceIndex, facej)] *
                             (phi[ijToid(index, j)] - phi[ijToid(index1, j)]) /
                             h[coord];
              if (i <= 1)
                values[0] = 1e4;
              if (i <= 0)
                values[1] = 1e4;
              if (i >= _gridSize[coord] - 2)
                values[5] = 1e4;
              if (i >= _gridSize[coord] - 1)
                values[4] = 1e4;
            } else {
              index = std::min(_gridSize[coord] - 1, std::max(1, j + inc + 1));
              index1 = std::min(_gridSize[coord] - 2, std::max(0, j + inc));
              faceIndex = std::min(_gridSize[coord], std::max(0, j + inc + 1));
              values[ival] = u[faceijToid(coord, facei, faceIndex)] *
                             (phi[ijToid(i, index)] - phi[ijToid(i, index1)]) /
                             h[coord];
              if (j <= 1)
                values[0] = 1e4;
              if (j <= 0)
                values[1] = 1e4;
              if (j >= _gridSize[coord] - 2)
                values[5] = 1e4;
              if (j >= _gridSize[coord] - 1)
                values[4] = 1e4;
            }
          }
        }

        double dPhi = Ramuh::Weno::evaluate(values, h[coord], false);
        weno[id] += dPhi;
      }
    }

    for (int id = 0; id < cellCount(); id++) {
      phi[id] = phi[id] - weno[id] * _dt;
    }

    _ellapsedTime += _dt;
  }

  void print() {
    auto &analytic = getScalarData(_analyticId);
    auto &phi = getScalarData(_functionId);
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

  double _dt, _ellapsedTime;
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
