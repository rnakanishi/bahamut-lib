#ifndef __RAMUH_HEAT_KERNEL_H__
#define __RAMUH_HEAT_KERNEL_H__

#include <Eigen/Sparse>
#include <map>
#include <geometry/bounding_box.h>

namespace Ramuh {

class VectorHeatKernel {
public:
  VectorHeatKernel(/* args */);

  VectorHeatKernel(Eigen::Array3d resolution);

  ~VectorHeatKernel();

  virtual void buildLaplacian();

  void addInitialCondition(int pid, Eigen::Array3d value);

  void resetInitialConditions();

  std::vector<Eigen::Array3d> solve();

protected:
  Eigen::Array3i _resolution;
  BoundingBox3 _domain;

  Eigen::SparseMatrix<double> L;
  Eigen::MatrixXd phi;

  std::map<int, Eigen::Array3d> _initialValues;

  double _dt;
};

} // namespace Ramuh

#endif