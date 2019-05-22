#ifndef __RAMUH_TRIANGLE_MESH_H__
#define __RAMUH_TRIANGLE_MESH_H__

#include <cmath>
#include <glm/vec3.hpp>
#include <map>
#include <vector>

namespace Ramuh {
class TriangleMesh {
public:
  TriangleMesh();
  ~TriangleMesh();

  /**
   * @brief Get the Vertices Size value
   *
   * @return uint
   **/
  uint getVerticesSize();

  /**
   * @brief Get the Faces Size
   *
   * @return uint
   **/
  uint getFacesSize();

  /**
   * @brief Add a vertex to the MeshModel structure and return its index
   *
   * @param vertex
   * @return uint
   **/
  virtual uint addVertex(glm::vec3 vertex);

  /**
   * @brief  Add a bundle of vertices and return their respective indices
   *
   * @return std::vector<uint>
   **/
  // virtual std::vector<uint> addVertices(std::vector<glm::vec3> vertices);

  /**
   * @brief Add A face, that corresponds to vertices indices that belong to face
   *
   * @param face
   * @return uint index of the added face
   **/
  virtual uint addFace(glm::ivec3 face);

  /**
   * @brief Add faces bundle, where each face has their corresponding vertices
   *
   * @param faces
   * @return std::vector<uint> indices of the added faces
   **/
  // virtual std::vector<uint> addFaces(glm::ivec3 faces);

  virtual void addVertexNormal(uint vertexId, glm::vec3 normal) {}

  // virtual void addFaceNormal(uint faceId, glm::vec3 normal);

  /**
   * @brief Remove very vertices with the same coordinates and adjust their
   * faces
   *
   * @param vertices
   * @return uint
   **/
  virtual uint removeDoubles() {}

  /**
   * @brief Get the vertex object corresponding to index
   *
   * @param index of the desired vertex
   * @return glm::vec3
   **/
  glm::vec3 getVertex(uint index);

  /**
   * @brief Get the Face object corresponding to index
   *
   * @param index of the desired face
   * @return glm::ivec3
   **/
  glm::ivec3 getFace(uint index);

  glm::vec3 getCentroid();

  glm::vec3 getBboxCenter();

  glm::vec3 getBBoxSize();

protected:
  void _computeCentroid();

  glm::vec3 _centroid;
  glm::vec3 _bboxMax, _bboxMin;

  std::vector<glm::vec3> _vertices;
  std::vector<glm::vec3> _vertexNormal;
  std::vector<glm::ivec3> _faces; //!< Face composed by three vertices
};

} // namespace Ramuh
#endif
