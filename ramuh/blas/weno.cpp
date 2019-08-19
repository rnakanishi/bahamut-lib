#include <blas/weno.h>
#include <Eigen/Dense>

namespace Ramuh {
Weno::Weno() {}

double Weno::evaluate(std::vector<double> values, double h, bool isNegative,
                      bool isBoundary) {
  Eigen::Vector3d w, phi, S, alpha;
  Eigen::VectorXd v(5); // values to form the convex combination

  if (!isNegative)
    for (int i = 0; i < 5; i++) {
      v[i] = (values[i + 1] - values[i]) / h;
    }
  else {
    for (int i = 0; i < 5; i++) {
      v[i] = (values[i] - values[i + 1]) / h;
    }
  }

  phi[0] = (2 * v[0] - 7 * v[1] + 11 * v[2]) / 6.;
  phi[1] = (-1 * v[1] + 5 * v[2] + 2 * v[3]) / 6.;
  phi[2] = (2 * v[2] + 5 * v[3] - 1 * v[4]) / 6.;

  double sConst = 13 / 12;
  S[0] = sConst * (v[0] - 2 * v[1] + v[2]) * (v[0] - 2 * v[1] + v[2]) +
         0.25 * (v[0] - 4 * v[1] + 3 * v[2]) * (v[0] - 4 * v[1] + 3 * v[2]);
  S[1] = sConst * (v[1] - 2 * v[2] + v[3]) * (v[1] - 2 * v[2] + v[3]) +
         0.25 * (v[1] - v[3]) * (v[1] - v[3]);
  S[2] = sConst * (v[2] - 2 * v[3] + v[4]) * (v[2] - 2 * v[3] + v[4]) +
         0.25 * (3 * v[2] - 4 * v[3] + v[4]) * (3 * v[2] - 4 * v[3] + v[4]);

  for (int i = 0; i < 5; i++) {
    if (std::fabs(v[i]) > 1000)
      v[i] = 0;
  }

  double eps =
      1e-6 * std::max(v[0] * v[0],
                      std::max(v[1] * v[1],
                               std::max(v[2] * v[2],
                                        std::max(v[3] * v[3], v[4] * v[4])))) +
      1e-32;
  alpha[0] = 0.1 / ((S[0] + eps) * (S[0] + eps));
  alpha[1] = 0.6 / ((S[1] + eps) * (S[1] + eps));
  alpha[2] = 0.3 / ((S[2] + eps) * (S[2] + eps));

  w[0] = alpha[0] / alpha.sum();
  w[1] = alpha[1] / alpha.sum();
  w[2] = alpha[2] / alpha.sum();

  return w.dot(phi);

  // TODO: implement smoothness function

} // namespace Ramuh

} // namespace Ramuh