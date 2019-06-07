#ifndef __GARUDA_LIGHT_HPP__
#define __GARUDA_LIGHT_HPP__

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <shader/shader.hpp>
#include <graphics/mesh_object.hpp>

namespace Garuda {

/**
 * @brief Basic Light class. Can also be usaed as a point light
 *
 */
class Light {
public:
  Light();

  glm::vec3 getPosition();

  glm::vec4 getColor();

  glm::mat4 &transformMatrix();

  void setPosition(glm::vec3 position);

  void setPosition(glm::vec4 color);

  void initialize();

  /**
   * @brief Draws a silhouette for the light source. As it is only a point, it
   * is represented by a small cube point
   *
   * @param shader lamp shader.
   */
  void draw(Shader shader);

  void illuminate(Shader shader);

  void move();

protected:
  glm::mat4 _transform;
  glm::vec3 _position;
  glm::vec4 _color;

  MeshObject _shape;
  unsigned int _vao;
};

/**
 * @brief Besides basic attributes from Light class, also includes a look-at
 * direction together with an opening angle for spot illumination
 *
 */
class SpotLight : public Light {
public:
  SpotLight();

  void setDirection(glm::vec3 direction);

  void setAngle(float angle);

  glm::vec3 getDirection();

  float getAngle();

protected:
  float _angle;
  glm::vec3 _direction;
};

} // namespace Garuda

#endif // !__GARUDA_LIGHT_HPP