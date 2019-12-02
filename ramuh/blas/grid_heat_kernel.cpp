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
  L.setZero();

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

std::vector<Eigen::Array3d>
GridHeatKernel::computeParallelTransport(int propagationDistance,
                                         bool normalize) {
  size_t nCells = _resolution[0] * _resolution[1] * _resolution[2];
  Eigen::SparseMatrix<double> M(nCells, nCells), Area(nCells, nCells);
  Eigen::Array3d h =
      _domain.getSize().cwiseQuotient(_resolution.cast<double>());
  double dt = h[0] * h[1] * h[2];
  Eigen::SparseMatrix<double> b(nCells, 3);
  Eigen::MatrixXd x;
  b.setFromTriplets(phiTriplets.begin(), phiTriplets.end());
  Area.setIdentity();
  Area *= h[0] * h[1];

  // for (size_t t = 0; t < propagationDistance; t++) {
  //   Eigen::SparseMatrix<double> aux = x;
  //   aux += L * x;
  //   x = aux;

  //   for (auto initial : _initialVectorValues) {
  //     x.coeffRef(initial.first, 0) = initial.second[0];
  //     x.coeffRef(initial.first, 1) = initial.second[1];
  //     x.coeffRef(initial.first, 2) = initial.second[2];
  //   }
  // }
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      cg;
  cg.compute(Area - L * dt);
  x = cg.solve(b);

  std::vector<Eigen::Array3d> result(nCells, Eigen::Array3d(0.0));

  // for (size_t out = 0; out < x.outerSize(); out++) {
  // for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(x,
  // out);
  //  it; ++it) {
  for (size_t it = 0; it < x.rows(); it++) {
    if (normalize)
      result[it] = x.row(it).normalized();
    else
      result[it] = x.row(it);
  }
  // }
  // }

  return result;
}

std::vector<double> GridHeatKernel::computeHeatDistance(
    std::vector<Eigen::Array3d> &gradientField) {

  int nCells = _resolution.prod();
  Eigen::Array3d h =
      _domain.getSize().cwiseQuotient(_resolution.cast<double>());

  if (gradientField.size() != nCells) {
    // TODO: Change to an exception
    std::cerr
        << "\033[1;31m[ERROR]\033[0m: "
        << " GridHeatKernel::computeHeatDistance: "
        << "gradientField has different size than the specified resolution\n";
    return std::vector<double>();
  }
  // Diffuse the heat from source points
  Eigen::SparseMatrix<double> Id(nCells, nCells);
  Eigen::SparseVector<double> bHeat(nCells);
  double dt = h[0] * h[0];
  Id.setIdentity();
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      cgGrad;
  cgGrad.compute(Id - L * dt);

  for (auto value : _initialScalarValues) {
    bHeat.coeffRef(value.first) = value.second;
  }
  Eigen::VectorXd xDistance = cgGrad.solve(bHeat);

  std::cerr << xDistance;
  // Normalize the gradient
  std::vector<double> distance;
  for (size_t cell = 0; cell < nCells; cell++) {
    distance.emplace_back(xDistance[cell]);
  }
  auto gradient = computeScalarFieldGradient(distance);

  // Vector field evaluation
  auto divergence = computeVectorFieldDivergence(gradient);
  Eigen::VectorXd b(nCells);
  for (size_t i = 0; i < nCells; i++) {
    b[i] = divergence[i];
  }
  // for (auto initCond : _initialScalarValues) {
  // b[initCond.first] = initCond.second;
  // L.row(initCond.first) *= 0; // Set entire row to 0
  // L.coeffRef(initCond.first, initCond.first) = 1;
  // }
  std::cerr << b << std::endl;

  Eigen::VectorXd x(nCells);

  // Laplacian linear system
  // Solving Poisson equation (Acutal distance computation)
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
  cg.compute(L / (h[0] * h[1]));
  x = cg.solve(b);
  std::cerr << L;

  std::vector<double> result(nCells, 0);

  for (size_t it = 0; it < x.size(); it++) {
    result[it] = x[it];
  }
  return result;
}

std::vector<double> GridHeatKernel::computeVectorFieldDivergence(
    std::vector<Eigen::Array3d> &gradient) {
  int nCells = _resolution.prod();
  std::vector<double> div(nCells, 0);
  Eigen::Array3d h =
      _domain.getSize().cwiseQuotient(_resolution.cast<double>());

  // Compute finite differences for all values
  for (size_t cellId = 0; cellId < nCells; cellId++) {
    auto ijk = idToijk(cellId);
    int i = ijk[0], j = ijk[1], k = ijk[2];

    for (size_t dim = 0; dim < 3; dim++) {
      double up = 0, down = 0, diff = 0;
      Eigen::Array3d direction(0);
      direction[dim] = 1;
      if (ijk[dim] > 0 && ijk[dim] < _resolution[dim] - 1) {
        Eigen::Array3d ijkArray(ijk[0], ijk[1], ijk[2]);
        Eigen::Array3d coords1 = ijkArray - direction;
        Eigen::Array3d coords2 = ijkArray + direction;
        diff = gradient[ijkToid(coords2[0], coords2[1], coords2[2])][dim] -
               gradient[ijkToid(coords1[0], coords1[1], coords1[2])][dim];
        diff /= 2 * h[dim];
      } else if (ijk[dim] > 0) {
        Eigen::Array3d ijkArray(ijk[0], ijk[1], ijk[2]);
        Eigen::Array3d coords = ijkArray - direction;
        diff = gradient[cellId][dim] -
               gradient[ijkToid(coords[0], coords[1], coords[2])][dim];
        diff /= h[dim];
      } else if (ijk[dim] < _resolution[dim] - 1) {
        Eigen::Array3d ijkArray(ijk[0], ijk[1], ijk[2]);
        Eigen::Array3d coords = ijkArray + direction;
        diff = gradient[ijkToid(coords[0], coords[1], coords[2])][dim] -
               gradient[cellId][dim];
        diff /= h[dim];
      }
      // Evaluate which direction to choose
      // if (ijk[dim] < _resolution[dim] - 1) {
      //   div[cellId] += up;
      // } else {
      //   div[cellId] += down;
      // }
      div[cellId] += diff;
    }
    // div[cellId] = gradient[cellId].sum() / 2;
  }

  return div;
}

std::vector<Eigen::Array3d>
GridHeatKernel::computeScalarFieldGradient(std::vector<double> &distance) {
  int nCells = _resolution.prod();
  std::vector<Eigen::Array3d> grad(nCells, Eigen::Array3d(0));
  Eigen::Array3d h =
      _domain.getSize().cwiseQuotient(_resolution.cast<double>());

  // Compute finite differences for all values
  for (size_t cellId = 0; cellId < nCells; cellId++) {
    auto ijk = idToijk(cellId);
    int i = ijk[0], j = ijk[1], k = ijk[2];

    for (size_t dim = 0; dim < 3; dim++) {
      double up = 0, down = 0, diff = 0;
      Eigen::Array3d direction(0);
      direction[dim] = 1;
      if (ijk[dim] > 0 && ijk[dim] < _resolution[dim] - 1) {
        Eigen::Array3d ijkArray(ijk[0], ijk[1], ijk[2]);
        Eigen::Array3d coords1 = ijkArray - direction;
        Eigen::Array3d coords2 = ijkArray + direction;
        diff = distance[ijkToid(coords2[0], coords2[1], coords2[2])] -
               distance[ijkToid(coords1[0], coords1[1], coords1[2])];
        diff /= 2 * h[dim];
      } else if (ijk[dim] > 0) {
        Eigen::Array3d ijkArray(ijk[0], ijk[1], ijk[2]);
        Eigen::Array3d coords = ijkArray - direction;
        diff = distance[cellId] -
               distance[ijkToid(coords[0], coords[1], coords[2])];
        diff /= h[dim];
      } else if (ijk[dim] < _resolution[dim] - 1) {
        Eigen::Array3d ijkArray(ijk[0], ijk[1], ijk[2]);
        Eigen::Array3d coords = ijkArray + direction;
        diff = distance[ijkToid(coords[0], coords[1], coords[2])] -
               distance[cellId];
        diff /= h[dim];
      }
      grad[cellId][dim] += diff;
    }
    grad[cellId] = grad[cellId].matrix().normalized().array();
    std::cerr << grad[cellId].transpose() << std::endl;
  }

  return grad;
}

} // namespace Ramuh