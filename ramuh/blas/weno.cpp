#include <blas/weno.h>
#include <Eigen/Dense>

namespace Ramuh {
Weno::Weno() {}

double Weno::evaluate(std::vector<double> values, double h, bool isLeft) {
  Eigen::Vector3d w, phi;
  Eigen::VectorXd v(5); // values to form the convex combination

  v[0] = (values[1] - values[0]) / h;
  v[1] = (values[2] - values[1]) / h;
  v[2] = (values[3] - values[2]) / h;
  v[3] = (values[4] - values[3]) / h;
  v[4] = (values[5] - values[4]) / h;

  phi[0] = (2 * v[0] - 7 * v[1] + 11 * v[2]) / 6.;
  phi[1] = (1 * v[1] + 5 * v[2] - 2 * v[3]) / 6.;
  phi[2] = (2 * v[2] + 5 * v[3] - 1 * v[4]) / 6.;

  w[0] = 0.1;
  w[1] = 0.6;
  w[2] = 0.3;

  return w.dot(phi);
}

} // namespace Ramuh