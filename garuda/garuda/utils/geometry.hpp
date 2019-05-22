#ifndef __GARUDA_UTILS_GEOMETRY_HPP
#define __GARUDA_UTILS_GEOMETRY_HPP
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

namespace Garuda {
class Geometry {
public:
  glm::mat4 axisAlign(glm::vec3 position, glm::vec3 direction);
}
} // namespace Garuda

#endif // !__GARUDA_UTILS_GEOMETRY_HPP
