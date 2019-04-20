#ifndef __GARUDA_TRIANGLE_MESH_HPP__
#define __GARUDA_TRIANGLE_MESH_HPP__

#include <glad/glad.h>
#include <structures/triangle_mesh.h>
#include <shader/shader.hpp>

namespace Garuda {
class MeshObject : public Ramuh::TriangleMesh {
public:
  MeshObject();

  /**
   * @brief Creates vbo, vao and ebo buffers.
   *
   **/
  void initialize();

  /**
   * @brief Bind corresponding buffers for rendering
   *
   **/
  void draw(Shader shader);

  /**
   * @brief Load an obj mesh and assemble buffers created before.
   *
   **/
  void loadObjMesh(const char *objPath);

  void addTextureCoordinate(glm::vec2 texCoord);

  uint addVertex(glm::vec3 vertex) override;

  void addVertexNormal(uint vertexId, glm::vec3 normal) override;

  glm::vec3 getCentroid();

  glm::vec3 getBboxCenter();

  glm::vec3 getBBoxSize();

protected:
  void _computeCentroid();

  /**
   * @brief This method do the reading function itself. Using tiny object
   *functions, assign properly the vertices coordinates and faces. If file
   *contains texture coordinates, they are assigned as well
   *
   * @param objPath
   **/
  void _readObj(const char *objPath);
  // TODO: change vbo to accept any size
  unsigned int _vbo[3], _vao, _ebo;
  std::vector<glm::vec3> _vertexColor;
  std::vector<glm::vec3> _vertexTexture;

  glm::vec3 _centroid;
  glm::vec3 _bboxMax, _bboxMin;
};

} // namespace Garuda

#endif