#include <geometry/geometry_utils.h>
#include <glm/geometric.hpp>

namespace Ramuh {

Geometry::Geometry() {}

glm::vec3 Geometry::closestPointTriangle(glm::vec3 p, glm::vec3 a, glm::vec3 b,
                                         glm::vec3 c) {
  glm::vec3 ab = b - a;
  glm::vec3 ac = c - a;
  glm::vec3 ap = p - a;
  float d1 = glm::dot(ab, ap);
  float d2 = glm::dot(ac, ap);

  // Check if P in edge AB
  glm::vec3 bp = p - b;
  float d3 = glm::dot(ab, bp);
  float d4 = glm::dot(ac, bp);
  float vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
    float v = d1 / (d1 - d3);
    return a + v * ab; // barycentric coordinates (1-v,v,0)
  }

  // Check if P in edge region of AC, if so return projection of P onto AC
  glm::vec3 cp = p - c;
  float d5 = glm::dot(ab, cp);
  float d6 = glm::dot(ac, cp);
  float vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    float w = d2 / (d2 - d6);
    return a + w * ac; // barycentric coordinates (1-w,0,w)
  }

  // Check if P in edge region of BC, if so return projection of P onto BC
  float va = d3 * d6 - d5 * d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return b + w * (c - b); // barycentric coordinates (0,1-w,w)
  }

  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  float denom = 1.0f / (va + vb + vc);
  float v = vb * denom;
  float w = vc * denom;
  return a + ab * v +
         ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
}

glm::vec3 Geometry::closestPointPlane(glm::vec3 p, glm::vec3 a, glm::vec3 b) {
  glm::vec3 ab = b - a;
  // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
  float t = glm::dot(p - a, ab) / glm::dot(ab, ab);
  // If outside segment, clamp t (and therefore d) to the closest endpoint
  if (t < 0.0f)
    t = 0.0f;
  if (t > 1.0f)
    t = 1.0f;
  // Compute projected position from the clamped t
  return a + t * ab;
}

} // namespace Ramuh