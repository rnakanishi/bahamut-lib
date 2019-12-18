#ifndef __RAMUH_STRUCTURES_LINE_MESH_HPP__
#define __RAMUH_STRUCTURES_LINE_MESH_HPP__

#include <Eigen/Dense>
#include <vector>
#include <map>

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
  int connectVertices(int vertex1, int vertex2);
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
   * @brief Get all the segments the vertex has.
   *
   * @param vertexId which segments are wnated
   * @return std::vector<int>& reference to a vector containing all vertex
   * connections
   */
  std::vector<int> &getVertexSegments(int vertexId);

  /**
   * @brief Get the reference for the Segments List, containing the pair of
   * vertices that are connected
   *
   * @return std::vector<Eigen::Array2i>&
   */
  std::vector<Eigen::Array2i> &getSegmentsList();

  /**
   * @brief Verify if the two given vertex ids are connected and return the id
   * of the segment connecting the vertices. If no connection is found, then -1
   * is returned instead
   *
   * @param vertex1
   * @param vertex2
   * @return the id for the segment if any. -1 if no connections was found
   */
  int hasConnection(int vertex1, int vertex2);

  /**
   * @brief Get the Number Of Connections for the given vertex Id. If no vertex
   * is found, or if no connections exist, then 0 is returned
   *
   * @param vertexId
   * @return int number of total connectios to the vertex
   */
  int getNumberOfConnections(int vertexId);

  std::vector<int> getAdjacentVertices(int vertexId);

  /**
   * @brief  Remove segment between both vertices. The segments list and each
   *vertex list are updated as well
   *
   * @param vertex1
   * @param vertex2
  **/
  void disconnectVertices(int vertex1, int vertex2);

protected:
  std::vector<Eigen::Array2d> _verticesPosition;
  std::vector<Eigen::Array2i> _segments; // Pair of connected vertices
  std::map<int, std::vector<int>>
      _vertSegments; // List of segments connecteds to the vertex
};

} // namespace Ramuh

#endif