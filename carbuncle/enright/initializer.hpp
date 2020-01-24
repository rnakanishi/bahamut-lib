#include <fluids/levelset_fluid2.h>
#include <fluids/levelset_fluid3.h>
#include <geometry/bounding_box.h>

namespace Carbuncle {
class LevelSetInitializer3 {

public:
  enum class Shape {
    DEFAULT = 0,
    SPHERE = 1,
    CUBE = 2,
    ZALESAK = 3,
    TORUS = 4
  };

  static void initializer(Leviathan::LevelSetFluid3 &levelset, Shape shape,
                          Eigen::Array3d center = Eigen::Array3d(0.0),
                          double radius = 1.0) {
    Ramuh::BoundingBox3 domain = levelset.getDomain();
    Eigen::Array3i resolution = levelset.getResolution();
    auto &phi = levelset.getCellScalarData("phi");

    switch (shape) {
    case Shape::SPHERE:
      for (size_t i = 0; i < levelset.cellCount(); i++) {
        phi[i] = std::sqrt(std::pow(p[0] - center[0], 2) +
                           std::pow(p[1] - center[1], 2)) -
                 radius;
      }
      break;

    case Shape::CUBE:
      for (size_t i = 0; i < levelset.cellCount(); i++) {
        Eigen::Array3d p = levelset.getCellPosition(i);
        Eigen::Array3d objCenter = p - center;
        double distance =
            std::max(std::max(std::fabs(objCenter[0]), std::fabs(objCenter[1])),
                     std::fabs(objCenter[2])) -
            radius;
        if (distance > 0) {
          p = p.abs();
          distance = 0.0;
          objCenter[0] = std::max(0.0, p[0] - radius);
          objCenter[1] = std::max(0.0, p[1] - radius);
          objCenter[2] = std::max(0.0, p[2] - radius);
          distance = objCenter.matrix().norm();
        }
        phi[i] = distance;
      }
      break;

    case Shape::ZALESAK:
      for (size_t i = 0; i < levelset.cellCount(); i++) {
        auto p = levelset.getCellPosition(i);
        if (pow(p[0] - center[0], 2) + pow(p[1] - center[1], 2) +
                pow(p[2] - center[2], 2) <
            radius * radius)
          phi[i] = 5;
        double s = radius / 4.;
        double s6 = radius / 6.;
        if (p[0] >= center[0] - s && p[0] <= center[0] + s &&
            p[1] < center[1] + 2 * s6 && p[1] > center[1] - radius - s6 &&
            p[2] >= center[2] - radius && p[2] <= center[2] + radius)
          phi[i] = -2;
      }
      break;

    case Shape::TORUS:
      for (size_t i = 0; i < levelset.cellCount(); i++) {
        Eigen::Array3d p = levelset.getCellPosition(i);
        p = p - center;
        p = p.cwiseProduct(p);
        double distance = (1.5 * radius - std::sqrt(p[0] + p[2])) *
                              (1.5 * radius - std::sqrt(p[0] + p[2])) +
                          p[1] - 0.75 * radius * radius;
      }
      break;

    default:
      break;
    }
  }
};

} // namespace Carbuncle