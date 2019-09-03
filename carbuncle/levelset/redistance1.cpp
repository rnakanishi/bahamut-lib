#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <structures/mac_grid1.h>

class RedistanceClass : public Ramuh::MacGrid1 {
public:
  RedistanceClass() : RedistanceClass(32, Ramuh::BoundingBox1::unitBox()) {}

  RedistanceClass(int gridSize, Ramuh::BoundingBox1 domain)
      : Ramuh::MacGrid1(domain, gridSize) {
    _phiId = newLabel("phi");
  }

  void initialize(double center, double radius) {
    auto &function = getScalarData(_phiId);
    auto h = getH();
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getPosition(i);

      function[i] = (p - 0.4 * h) * (p + 6) / 2 + 1;
    }
  }

  void defineVelocity() {}

  void print() {
    auto &phi = getScalarData(_phiId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/redistance/1d/" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);

    for (size_t i = 0; i < cellCount(); i++) {
      file << phi[i] << " ";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

  void redistance() {
    auto h = getH();
    auto &phi = getScalarData(_phiId);
    double eps = h, error = 1e8, dt = 0.5 * h;
    std::vector<double> gradient(cellCount());
    std::vector<double> cellSignal(cellCount());
    std::vector<double> initialPhi(cellCount());
    std::vector<double> interfaceFactor(cellCount());
    std::vector<bool> isInterface(cellCount(), false);

    for (int id = 0; id < cellCount(); id++) {
      initialPhi[id] = phi[id];
      if (phi[id] > 0)
        cellSignal[id] = 1;
      else if (phi[id] < 0)
        cellSignal[id] = -1;
      else
        cellSignal[id] = 0;

      isInterface[id] = false;
      interfaceFactor[id] = 0;
      if (id > 0 && id < _gridSize - 1) {
        if (phi[id] * phi[id - 1] <= 0 || phi[id] * phi[id + 1] <= 0) {
          isInterface[id] = true;
          interfaceFactor[id] = h * 2. * phi[id];
          interfaceFactor[id] /= abs(phi[id + 1] - phi[id - 1]);
        }
        // cellSignal[id] = interfaceFactor[id] / h;
      }
    }

    int t;
    for (t = 0; t < 1; t++) {
      // #pragma omp parallel for
      for (int id = 0; id < cellCount(); id++) {
        int i = id;
        double dx[2]; // Index 0: Dx-, Index 1: Dx+
        dx[0] = dx[1] = 0;

        if (i > 0)
          dx[0] = (phi[id] - phi[id - 1]) / h;
        if (i < _gridSize - 1)
          dx[1] = (phi[id + 1] - phi[id]) / h;

        double gx;
        if (initialPhi[id] > 0) {
          double a, b;
          a = std::max(0., dx[0]);
          b = std::min(0., dx[1]);
          gx = std::max(abs(a), abs(b));
        } else if (initialPhi[id] < 0) {
          double a, b;
          a = std::min(0., dx[0]);
          b = std::max(0., dx[1]);
          gx = std::max(abs(a), abs(b));
        } else {
          gx = 1;
        }
        gradient[id] = abs(gx) - 1;
      }

      error = 0.;
      // #pragma omp parallel for reduction(+ : error)
      for (int id = 0; id < cellCount(); id++) {
        double newPhi = 0.;

        if (!isInterface[id]) {
          // if (true) {
          newPhi = phi[id] - dt * cellSignal[id] * gradient[id];
        } else {
          newPhi = phi[id] - (dt / h) * (cellSignal[id] * abs(phi[id]) -
                                         interfaceFactor[id]);
        }
        error += abs(phi[id] - newPhi);

        phi[id] = newPhi;
      }
      if (error / cellCount() < dt * h * h)
        break;
    }
  }

protected:
  int _phiId;
};

int main(int argc, char const *argv[]) {
  RedistanceClass cubes(32, Ramuh::BoundingBox1(-5, 5));

  cubes.initialize(0.0, 4);

  cubes.print();

  for (int i = 1; i <= 300; i++) {
    cubes.redistance();

    cubes.print();
  }

  //   auto surface = cubes.marc hingTetrahedra();
  //   Ramuh::FileWriter writer;
  //   std::ostringstream objname;
  //   objname << "results/marching/tetra_" << std::setfill('0') << std::setw(4)
  //   << 0
  //   << ".obj";
  //   writer.writeMeshModel(surface, objname.str());
  // }
  return 0;
}
