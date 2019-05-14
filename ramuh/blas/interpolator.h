#ifndef __RAMUH_BLAS_INTERPOLATOR_H
#define __RAMUH_BLAS_INTERPOLATOR_H

#include <vector>

namespace Ramuh {
class Interpolator {
public:
  Interpolator();

  static double linear(double position, std::vector<double> points,
                       std::vector<double> values);
};
} // namespace Ramuh

#endif