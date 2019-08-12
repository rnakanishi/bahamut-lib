#ifndef __RAMUH_BLAS_WENO_H__
#define __RAMUH_BLAS_WENO_H__
#include <vector>

namespace Ramuh {

class Weno {
public:
  Weno();

  /**
   * @brief Uses values given as input to evaluate diferential values at central
   * point. For left differencing, vector index [3] indicate such central point,
   * otherwise, if right differencing is used, then vector indes [2] is the
   * central point. The third parameter controls if the upwind velocity is
   * coming from left or right: if true, computes left differencig.
   *
   * Reference: Stanley Osher, Ronald Fedkiw. "Level set methods and dynamic
   * implicit surfaces", Springer 2003
   *
   * @param values points over values to be computed
   * @param h grid spacing to be used in the computations
   * @param isLeft controls which side should be use on differentiation
   * @return double differential value computed
   */
  static double evaluate(std::vector<double> values, double h,
                         bool isLeft = true);

private:
};

} // namespace Ramuh

#endif