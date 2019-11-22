#include <blas/heat_kernel.h>

namespace Ramuh {

VectorHeatKernel::VectorHeatKernel(/* args */) {}

VectorHeatKernel::~VectorHeatKernel() {}

std::vector<Eigen::Array3d> VectorHeatKernel::solve() {}

void VectorHeatKernel::buildLaplacian() {}

void VectorHeatKernel::addInitialCondition(int id, Eigen::Array3d value) {

  // L.row(id) *= 0; // Set entire row to 0
  // L.coeffRef(id, id) = 1;

  phi.row(id) << value[0], value[1], value[2];
}

} // namespace Ramuh
