#include <shader/material.hpp>

namespace Garuda {
Material::Material() {}

Material::Material(glm::vec4 ka, glm::vec4 kd, glm::vec4 ks, float shininess) {
  _ka = ka;
  _kd = kd;
  _ks = ks;
  _shininess = shininess;
}

glm::vec4 Material::getAmbient() { return _ka; }

glm::vec4 Material::getDiffuse() { return _kd; }

glm::vec4 Material::getSpecular() { return _ks; }

float Material::getShininess() { return _shininess; }
} // namespace Garuda