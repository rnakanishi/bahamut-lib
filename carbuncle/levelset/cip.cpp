#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <structures/mac_grid1.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

class RedistanceClass : public Ramuh::MacGrid1 {
public:
  RedistanceClass() : RedistanceClass(32, Ramuh::BoundingBox1::unitBox()) {}

  RedistanceClass(int gridSize, Ramuh::BoundingBox1 domain)
      : Ramuh::MacGrid1(domain, gridSize) {
    _phiId = newLabel("phi");
    _gradientId = newLabel("gradient");
    _velocityId = newLabel("velocity");
  }

  void initialize(double center, double radius) {
    auto &phi = getScalarData(_phiId);
    auto &gradient = getScalarData(_gradientId);
    auto h = getH();
    for (size_t i = 0; i < cellCount(); i++) {
      auto p = getPosition(i);
      phi[i] = 0.;
      if (p >= -4. && p <= -3.)
        phi[i] = 5;
    }
    for (size_t i = 1; i < cellCount(); i++) {
      gradient[i] = (phi[i] - phi[i - 1]) / h;
    }
  }

  void defineVelocity() {
    auto &velocity = getScalarData(_velocityId);
    for (size_t i = 0; i < cellCount(); i++) {
      velocity[i] = 1.;
    }
  }

  void print() {
    auto &phi = getScalarData(_phiId);
    static int count = 0;
    std::ofstream file;
    std::stringstream filename;
    filename << "results/cip/1d/" << count++;
    file.open(filename.str().c_str(), std::ofstream::out);

    for (size_t i = 0; i < cellCount(); i++) {
      file << phi[i] << " ";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

  void advectionCip() {
    auto h = getH();
    auto &phi = getScalarData(_phiId);
    auto &gradient = getScalarData(_gradientId);
    auto &velocity = getScalarData(_velocityId);
    double eps = h, error = 1e8, dt = 0.5 * h;
    double cdt, kappa, kappa2;
    static int count = 0;
    std::vector<double> a(cellCount()), b(cellCount());
    std::vector<double> newGrad(cellCount()), integral(cellCount()),
        Dft((cellCount()));

    cdt = -velocity[0] * dt;
    kappa = -cdt / h;
    kappa2 = kappa * kappa;

    newGrad[0] = 0;
    Eigen::VectorXd newPhi;
    {
      for (size_t i = 0; i < cellCount() - 1; i++) {
        // Integral of f and Dfi
        integral[i] = phi[i]; //(phi[i - 1] + 2 * phi[i] + phi[i + 1]) / 4 * h;

        Dft[i] = (-kappa / 8 + kappa2 / 8 + kappa * kappa2 / 6 -
                  kappa2 * kappa2 / 4) *
                 gradient[i + 1] * h * h;
        Dft[i] += (kappa / 8 + kappa2 / 8 - kappa2 * kappa2 / 6 -
                   kappa2 * kappa2 / 2) *
                  gradient[i] * h * h;
        Dft[i] += (0.5 * kappa - 0.75 * kappa2 + 0.5 * kappa2 * kappa2) *
                  phi[i + 1] * h;
        Dft[i] +=
            (0.5 * kappa + 0.75 * kappa2 - 0.5 * kappa2 * kappa2) * phi[i] * h;

        // Predict next t gradient: Shifted gradients
        // Interpolate gradient  f'(xi, tn+1) = f'(xi - c*dt, tn);
        newGrad[i] = (gradient[i] - gradient[i - 1]) * (-cdt) /
                     (getPosition(i) - getPosition(i - 1));
      }
      // Calculate fn+1 (Equation 6a)
      Eigen::SparseMatrix<double> A(cellCount(), cellCount());
      Eigen::VectorXd b(cellCount());
      std::vector<Eigen::Triplet<double>> triplets;
      b[0] = b[cellCount() - 1] = 0;
      triplets.emplace_back(0, 0, 1);
      triplets.emplace_back(cellCount() - 1, cellCount() - 1, 1);
      for (size_t i = 1; i < cellCount() - 1; i++) {
        triplets.emplace_back(i, i + 1, 18. / 192.);
        triplets.emplace_back(i, i + 0, 156 / 192.);
        triplets.emplace_back(i, i - 1, 18. / 192);
        b[i] = (5. / 192) * (newGrad[i + 1] - newGrad[i - 1]) * h;
        b[i] += (1. / 192) * (18 * integral[i + 1] + 156 * integral[i] +
                              18 * integral[i - 1] - 5 * gradient[i + 1] * h +
                              5 * gradient[i - 1] * h);
        b[i] -= (Dft[i] - Dft[i - 1]) / h;
      }
      A.setFromTriplets(triplets.begin(), triplets.end());
      Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
      solver.analyzePattern(A);
      solver.factorize(A);
      // std::cerr << A;
      // std::cerr << b;

      // solver.compute(A);
      if (solver.info() != Eigen::Success)
        std::cerr << "Solver error: compute(A) failed\n";
      newPhi = solver.solve(b);
      if (solver.info() != Eigen::Success)
        std::cerr << "Solver error: solve(b) failed\n";

      count++;
      if (count < 10 || count % 50 == 0)
        // Recalculate f'n+1 from fn+1 using minor gradients
        for (size_t i = 1; i < cellCount() - 1; i++) {
          double left, right;
          left = (newPhi[i] - newPhi[i - 1]) / h;
          right = (newPhi[i + 1] - newPhi[i]) / h;
          if (abs(left) < abs(right))
            newGrad[i] = left;
          else
            newGrad[i] = right;
        }
    }
    for (size_t i = 1; i < cellCount() - 1; i++) {
      phi[i] = newPhi[i];
      gradient[i] = newGrad[i];
    }
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
  int _gradientId;
  int _velocityId;
};

int main(int argc, char const *argv[]) {
  RedistanceClass cubes(1000, Ramuh::BoundingBox1(-5, 5));

  cubes.initialize(0.0, 4);
  cubes.defineVelocity();

  cubes.print();

  for (int i = 1; i <= 1000; i++) {
    cubes.advectionCip();

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
