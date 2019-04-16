#ifndef __RAMUH_MESH_MODEL_H__
#define __RAMUH_MESH_MODEL_H__

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <glm/vec3.hpp>

namespace Ramuh {
class MeshModel3 {
public:
  MeshModel3();
  ~MeshModel3();

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
  uint addVertex(glm::vec3 vertex);

  /**
   * @brief  Add a bundle of vertices and return their respective indices
   *
   * @return std::vector<uint>
   **/
  std::vector<uint> addVertices(std::vector<glm::vec3> vertices);

  /**
   * @brief Add A face, that corresponds to vertices indices that belong to face
   *
   * @param face
   * @return uint index of the added face
   **/
  uint addFace(glm::ivec3 face);

  /**
   * @brief Add faces bundle, where each face has their corresponding vertices
   *
   * @param faces
   * @return std::vector<uint> indices of the added faces
   **/
  std::vector<uint> addFaces(glm::ivec3 faces);

  /**
   * @brief Remove very vertices with the same coordinates and adjust their
   * faces
   *
   * @param vertices
   * @return uint
   **/
  uint removeDoubles();

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

protected:
  class vec3Compare {
  public:
    bool operator()(const glm::vec3 v1, const glm::vec3 v2) const {
      return v1[0] < v2[0] ||
             (std::fabs(v1[0] - v2[0]) < 1e-8 && v1[1] < v2[1]) ||
             (std::fabs(v1[0] - v2[0]) < 1e-8 &&
              std::fabs(v1[1] - v2[1]) < 1e-8 && v1[2] < v2[2]);
    }
  };
  std::map<glm::vec3, uint, vec3Compare> _vMap;
  std::vector<glm::vec3> _vertices;
  std::vector<glm::ivec3> _faces; //!< Face composed by three vertices
  std::vector<std::vector<uint>> _usedVertices; //!< Map which faces use which
                                                //!< vertices. First index
                                                //!< correspond to vertex.
};

} // namespace Ramuh
#endif
