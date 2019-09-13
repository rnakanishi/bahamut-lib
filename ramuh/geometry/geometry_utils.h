#ifndef __RAMUH_GEOMETRY_UTILS_H__
#define __RAMUH_GEOMETRY_UTILS_H__
#include <glm/vec3.hpp>
#include <glm/vec2.hpp>

namespace Ramuh {

class Geometry {
public:
  Geometry();

  static glm::vec3 closestPointTriangle(glm::vec3 P, glm::vec3 a, glm::vec3 b,
                                        glm::vec3 c);

  static glm::vec3 closestPointPlane(glm::vec3 p, glm::vec3 a, glm::vec3 b);
  static glm::vec2 closestPointPlane(glm::vec2 p, glm::vec2 a, glm::vec2 b);
};

} // namespace Ramuh

#endif