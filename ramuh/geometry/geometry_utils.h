#ifndef __RAMUH_GEOMETRY_UTILS_H__
#define __RAMUH_GEOMETRY_UTILS_H__
#include <glm/vec3.hpp>

namespace Ramuh {

class Geometry {
public:
  Geometry();

  static glm::vec3 closestPointTriangle(glm::vec3 P, glm::vec3 a, glm::vec3 b,
                                        glm::vec3 c);

  static glm::vec3 closestPointPlane(glm::vec3 p, glm::vec3 a, glm::vec3 b);
};

} // namespace Ramuh

#endif