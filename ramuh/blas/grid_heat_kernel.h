#ifndef __RAMUH_GRID_HEAT_KERNEL_H__
#define __RAMUH_GRID_HEAT_KERNEL_H__

#include <blas/heat_kernel.h>
#include <geometry/bounding_box.h>

namespace Ramuh {
class GridHeatKernel : public VectorHeatKernel {
public:
  GridHeatKernel(/* args */);

  GridHeatKernel(Eigen::Array3i resolution, BoundingBox3 domain);

  std::vector<size_t> idToijk(size_t id);

  size_t ijkToid(size_t i, size_t j, size_t k);

  void buildLaplacian() override;

  std::vector<Eigen::Array3d> solve();

  ~GridHeatKernel();

protected:
};

} // namespace Ramuh

#endif