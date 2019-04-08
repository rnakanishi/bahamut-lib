#ifndef __RAMUH_MESH_MODEL_H__
#define __RAMUH_MESH_MODEL_H__

#include <Eigen/Dense>
#include <vector>

namespace Ramuh {
class MeshModel3 {
public:
  MeshModel3();
  ~MeshModel3();

  /**
   * @brief Add a vertex to the MeshModel structure and return its index
   *
   * @param vertex
   * @return uint
   **/
  uint addVertex(Eigen::Vector3d vertex);

  /**
   * @brief  Add a bundle of vertices and return their respective indices
   *
   * @return std::vector<uint>
   **/
  std::vector<uint> addVertices(std::vector<Eigen::Vector3d> vertices);

  /**
   * @brief Add A face, that corresponds to vertices indices that belong to face
   *
   * @param face
   * @return uint index of the added face
   **/
  uint addFace(Eigen::Vector3i face);

  /**
   * @brief Add faces bundle, where each face has their corresponding vertices
   *
   * @param faces
   * @return std::vector<uint> indices of the added faces
   **/
  std::vector<uint> addFaces(Eigen::Vector3i faces);

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
   * @return Eigen::Vector3d
   **/
  Eigen::Vector3d getVertex(uint index);

  /**
   * @brief Get the Face object corresponding to index
   *
   * @param index of the desired face
   * @return Eigen::Vector3i
   **/
  Eigen::Vector3i getFace(uint index);

protected:
  std::vector<Eigen::Vector3d> _vertices;
  std::vector<Eigen::Vector3i> _faces; //!< Face composed by three vertices
  std::vector<std::vector<uint>> _usedVertices; //!< Map which faces use which
                                                //!< vertices. First index
                                                //!< correspond to vertex.
};

} // namespace Ramuh
#endif
