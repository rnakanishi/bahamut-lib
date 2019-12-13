#ifndef __RAMUH_STRUCTURES_LINE_MESH_HPP__
#define __RAMUH_STRUCTURES_LINE_MESH_HPP__

#include <Eigen/Dense>
#include <vector>

namespace Ramuh {

class LineMesh {
public:
  LineMesh();

  /**
   * @brief Add a vertex to the mesh. Multiple vertices can be added as well
   *
   * @param position of the vertex to be added. if a std::vector is passed
   * instead, then multiple vertices are added to the structure
   *
   * @return return the indices of added vertices
   */
  int addVertex(Eigen::Array2d position);
  std::vector<int> addVertices(std::vector<Eigen::Array2d> positions);

  /**
   * @brief Add two vertices: one at @origin and the other at @ending position.
   * A connection is made between those vertices, adding a segment to the
   * structure, also.
   *
   * Vertices are expected to be given in CCW direction, otherwise, it may
   * affect normal generation
   *
   * @param origin Starting vertex of the segment
   * @param ending Ending vertex of the segment
   *
   * @return the index of the added segment
   */
  int addSegment(Eigen::Array2d origin, Eigen::Array2d ending);

  /**
   * @brief Connect two vertices id to form a segment. If a vector containing
   * the pair of vertices are passed as parameter, then multiple segments are
   * created
   *
   * @param segment vertices id to be connected. A vector containing the vertex
   * ids can be used as well
   *
   * @return the new segment id
   */
  int connectVertices(Eigen::Array2i segment);
  std::vector<int> connectVertices(std::vector<Eigen::Array2i> segment);

  /**
   * @brief Get the Vertex Position for the given id.
   *
   * @param VertexId of the vertex
   * @return Eigen::Array2d position in the space of the vertex
   */
  Eigen::Array2d getVertexPosition(int VertexId);

  /**
   * @brief Get the Vertices ids belonging to the segment. Vertices are returned
   * in CCW direction (order they are connected)
   *
   * @param segmentId
   * @return Eigen::Array2i vertices id
   */
  Eigen::Array2i getSegmentVertices(int segmentId);

  /**
   * @brief Get the Vertices count
   *
   * @return int number of vertices in the structure
   */
  int getVerticesCount();

  /**
   * @brief Get the Segments count
   *
   * @return int number of segments in the structure
   */
  int getSegmentsCount();

  /**
   * @brief Get the reference for the Vertices List, containing all the vertices
   * position
   *
   * @return std::vector<Eigen::Array2d>&
   */
  std::vector<Eigen::Array2d> &getVerticesList();

  /**
   * @brief Get the reference for the Segments List, containing the pair of
   * vertices that are connected
   *
   * @return std::vector<Eigen::Array2i>&
   */
  std::vector<Eigen::Array2i> &getSegmentsList();

protected:
  std::vector<Eigen::Array2d> _verticesPosition;
  std::vector<Eigen::Array2i> _segments;
};

} // namespace Ramuh

#endif