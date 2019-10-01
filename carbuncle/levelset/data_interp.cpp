#include <fluids/levelset_fluid2.h>
#include <cstdlib>
#include <iostream>

class InterpolationTest : public Leviathan::LevelSetFluid2 {
public:
  InterpolationTest()
      : Leviathan::LevelSetFluid2(Eigen::Array2i(64),
                                  Ramuh::BoundingBox2(-5, 5)) {}

  void initializeData() {
    auto &phi = getCellScalarData(_phiId);
    for (size_t i = 0; i < cellCount(); i++) {
      Eigen::Array2d pos = getCellPosition(i);
      phi[i] = sin(3 * pos[0]) * cos(3 * pos[1]);
    }
  }

  void interpolate() {
    auto h = getH();
    for (size_t i = 0; i < 1000; i++) {
      Eigen::Array2d target;
      Eigen::Array2d randFactor((std::rand() % 10000000) / 10000000.,
                                (std::rand() % 10000000) / 10000000.);
      target = _domain.min() +
               ((_domain.max() - _domain.min()).cwiseProduct(randFactor));
      double data = interpolateCellScalarData(_phiId, target);
      double analytic = sin(3 * target[0]) * cos(3 * target[1]);
      double error = (data - analytic);
      double order = log(std::abs(error)) / log(h[0]);

      std::cerr << target[0] << " " << target[1] << "\t";
      std::cerr << data << " " << analytic << "\t";
      std::cerr << error << " " << order << std::endl;
    }
  }

protected:
};

int main(int argc, char const *argv[]) {
  InterpolationTest cube;

  cube.initializeData();
  cube.interpolate();

  return 0;
}
