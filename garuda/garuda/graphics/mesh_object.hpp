#ifndef __GARUDA_TRIANGLE_MESH_HPP__
#define __GARUDA_TRIANGLE_MESH_HPP__

#include <glad/glad.h>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <map>
#include <shader/shader.hpp>
#include <shader/material.hpp>
#include <shader/texture.hpp>
#include <structures/triangle_mesh.h>
#include <graphics/scene_object.hpp>

namespace Garuda {
class MeshObject : public SceneObject, public Ramuh::TriangleMesh {
public:
  MeshObject();

  /**
   * @brief Creates vbo, vao and ebo buffers.
   *
   **/
  void initialize() override;

  /**
   * @brief Bind corresponding buffers for rendering
   *
   **/
  virtual void draw(Shader shader) override;

  /**
   * @brief Load an obj mesh and assemble buffers created before.
   *
   **/
  void loadObjMesh(const char *objPath);

  void centerizeObject();

  void loadTexture();

  void addTextureCoordinate(glm::vec2 texCoord);

  uint addVertex(glm::vec3 vertex) override;

  void assignVertices(std::vector<glm::vec3> &vertices);

  void assignNormals(std::vector<glm::vec3> &normals);

  void removeAllVertices();

  void addVertexNormal(uint vertexId, glm::vec3 normal) override;

  Texture &getTexture();

  bool &hasTexture();
  bool &hasNormal();
  bool &hasMaterial();

  void hasTexture(bool value);
  void hasNormal(bool value);
  void hasMaterial(bool value);

protected:
  /**
   * @brief This method do the reading function itself. Using tiny object
   *functions, assign properly the vertices coordinates and faces. If file
   *contains texture coordinates, they are assigned as well
   *
   * @param objPath
   **/
  void _readObj(const char *objPath);
  // TODO: change vbo to accept any size
  std::vector<glm::vec3> _vertexColor;
  std::vector<glm::vec2> _vertexTexture;

  Texture _textureImage;
  Material _material;
  bool _hasTexture, _hasNormal, _hasMaterial;
};

} // namespace Garuda

#endif