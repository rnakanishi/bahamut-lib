#include <blas/grid_heat_kernel.h>
#include <omp.h>
#include <iostream>
namespace Ramuh {

GridHeatKernel::GridHeatKernel(/* args */) {}

GridHeatKernel::GridHeatKernel(Eigen::Array3i resolution, BoundingBox3 domain) {
  _resolution = resolution;
  _domain = domain;
  this->buildLaplacian();
}
GridHeatKernel::~GridHeatKernel() {}

std::vector<size_t> GridHeatKernel::idToijk(size_t id) {
  std::vector<size_t> indices(3);
  indices[2] = id / (_resolution[0] * _resolution[1]);
  indices[1] = (id % (_resolution[0] * _resolution[1])) / _resolution[0];
  indices[0] = (id % (_resolution[0] * _resolution[1])) % _resolution[0];
  return indices;
}

size_t GridHeatKernel::ijkToid(size_t i, size_t j, size_t k) {
  return k * _resolution[0] * _resolution[1] + j * _resolution[0] + i;
}

void GridHeatKernel::buildLaplacian() {
  size_t nCells = _resolution[0] * _resolution[1] * _resolution[2];
  Eigen::Array3d h =
      _domain.getSize().cwiseQuotient(_resolution.cast<double>());
  std::vector<Eigen::Triplet<double>> triplets;

  L = Eigen::SparseMatrix<double>(nCells, nCells);
  phi = Eigen::MatrixXd(nCells, 3);
  L.setZero();
  phi.setZero();

#pragma omp parallel
  {
    std::vector<Eigen::Triplet<double>> threadTriplet;
    threadTriplet.clear();
#pragma omp for
    for (size_t id = 0; id < nCells; id++) {
      auto ijk = idToijk(id);
      size_t i = ijk[0], j = ijk[1], k = ijk[2];

      if (i > 0) {
        threadTriplet.emplace_back(id, ijkToid(i - 1, j, k), 1 / (h[0] * h[0]));
        threadTriplet.emplace_back(id, id, -1 / (h[0] * h[0]));
      }
      if (i < _resolution[0] - 1) {
        threadTriplet.emplace_back(id, ijkToid(i + 1, j, k), 1 / (h[0] * h[0]));
        threadTriplet.emplace_back(id, id, -1 / (h[0] * h[0]));
      }
      if (j > 0) {
        threadTriplet.emplace_back(id, ijkToid(i, j - 1, k), 1 / (h[1] * h[1]));
        threadTriplet.emplace_back(id, id, -1 / (h[1] * h[1]));
      }
      if (j < _resolution[1] - 1) {
        threadTriplet.emplace_back(id, ijkToid(i, j + 1, k), 1 / (h[1] * h[1]));
        threadTriplet.emplace_back(id, id, -1 / (h[1] * h[1]));
      }
      if (k > 0) {
        threadTriplet.emplace_back(id, ijkToid(i, j, k - 1), 1 / (h[2] * h[2]));
        threadTriplet.emplace_back(id, id, -1 / (h[2] * h[2]));
      }
      if (k < _resolution[2] - 1) {
        threadTriplet.emplace_back(id, ijkToid(i, j, k + 1), 1 / (h[2] * h[2]));
        threadTriplet.emplace_back(id, id, -1 / (h[2] * h[2]));
      }
    }
#pragma omp critical
    {
      triplets.insert(triplets.end(), threadTriplet.begin(),
                      threadTriplet.end());
    }
  }
  L.setFromTriplets(triplets.begin(), triplets.end());
}

std::vector<Eigen::Array3d> GridHeatKernel::solve() {
  size_t nCells = _resolution[0] * _resolution[1] * _resolution[2];
  Eigen::SparseMatrix<double> M(nCells, nCells);
  Eigen::Array3d h =
      _domain.getSize().cwiseQuotient(_resolution.cast<double>());
  double dt = h[0] * h[1] * h[2];
  Eigen::MatrixXd x;
  x = phi;

  for (size_t t = 0; t < 10; t++) {
    x = x + (dt * L) * x;
  }

  std::vector<Eigen::Array3d> result;
  for (size_t id = 0; id < nCells; id++) {
    result.emplace_back(
        Eigen::Array3d(x.row(id)[0], x.row(id)[1], x.row(id)[2]));
  }
  return result;
}

} // namespace Ramuh