#ifndef __GARUDA_MATERIAL_HPP__
#define __GARUDA_MATERIAL_HPP__
#include <glm/vec4.hpp>

namespace Garuda {
class Material {
public:
  Material();

  Material(glm::vec4 ka, glm::vec4 kd, glm::vec4 ks, float shininess);

  glm::vec4 getAmbient();

  glm::vec4 getDiffuse();

  glm::vec4 getSpecular();

  float getShininess();

private:
  glm::vec4 _ka, _kd, _ks;
  float _shininess;
};

} // namespace Garuda

#endif
