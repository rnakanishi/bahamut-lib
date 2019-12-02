#include <blas/heat_kernel.h>

namespace Ramuh {

VectorHeatKernel::VectorHeatKernel(/* args */) {}

VectorHeatKernel::~VectorHeatKernel() {}

std::vector<Eigen::Array3d> VectorHeatKernel::computeParallelTransport() {
  return std::vector<Eigen::Array3d>();
}

void VectorHeatKernel::buildLaplacian() {}

void VectorHeatKernel::addInitialCondition(int id, double value) {
  _initialScalarValues[id] = value;
}

void VectorHeatKernel::addInitialCondition(int id, Eigen::Array3d value) {

  // L.row(id) *= 0; // Set entire row to 0
  // L.coeffRef(id, id) = 1;

  // phi.row(id) << value[0], value[1], value[2];
  for (int i = 0; i < 3; i++) {
    phiTriplets.emplace_back(id, i, value[i]);
  }

  _initialVectorValues[id] = value;
}

} // namespace Ramuh
