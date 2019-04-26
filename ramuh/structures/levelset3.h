#ifndef __RAMUH_LEVELSET3_H__
#define __RAMUH_LEVELSET3_H__

#include <structures/grid3.h>
#include <structures/triangle_mesh.h>
#include <utils/material.h>
#include <glm/glm.hpp>

namespace Ramuh {

class LevelSet3 : public RegularGrid3 {
public:
  LevelSet3();

  ///
  /// Add an implicit region for the fluid
  void addImplicitFunction();

  ///
  /// Level set often gets unusual beahvior when advected. To avoid unphysical
  /// behavior, a redistancing function is applied so the level set keep its
  /// property as signed distance
  void redistance();

  ///
  /// Advect level set according to \f$\phi_t + u\cdot\nabla\phi = 0\f$. This
  /// method assumes that a velocity field is already defined over cell corners
  /// TODO: Change to semi lagrangean method
  void integrateLevelSet();

  /**
   * @brief Perform Mac Cormack advection method over levelset values
   *
   **/
  void macCormackAdvection();

  void solvePressure() override;

  ///
  /// Define a value for each vertex of the grid correspnoding to the isocontour
  /// of a sphere given its center and radius
  /// \param center center of the sphere
  /// \param radius radius of the sphere
  void addSphereSurface(Vector3d center, double radius);
  void addCubeSurface(Vector3d lower, Vector3d upper);

  /**
   * @brief Perform marching tetrahedra algorithm looking for zero level set.
   * The coordinates of the vertices correspond to the cell center. Possibly,
   * created vertices may be doubled and faces completely independent
   *
   * @return TriangleMesh full model with vertices and faces of the zero
   *levelset
   **/
  TriangleMesh marchingTetrahedra();

  void setPressureSecondOrder(bool value);

  void checkCellMaterial();

  void printVertexVelocity();

  void printLevelSetValue();

  void setResolution(Vector3i newResolution) override;

  std::vector<std::vector<double>> &operator[](const int i);

  Vector3d operator()(int i, int j, int k);

protected:
  /**
   * @brief Solve Eikonal equation for a given cell index. Look for all
   * direction neighbors choosing those that are nearest to the surface (in each
   * direction). solve linear, triangular or tetrahedral distance with chosen
   * neighbors
   *
   * @param cell cell integer index
   * @return double distance value to the levelset
   **/
  double _solveEikonal(glm::ivec3 cell);

  /**
   * @brief Auxiliary functino for marching tetrahedra
   * TODO: Check parameters and return value accordingly
   **/
  void _triangulate(std::vector<glm::ivec3> vertices, TriangleMesh &mesh);

  double _interpolatePhi(Eigen::Array3d position);

  glm::vec3 _findSurfaceCoordinate(glm::ivec3 v1, glm::ivec3 v2);
  Eigen::Array3d _findSurfaceCoordinate(Eigen::Array3i v1, Eigen::Array3i v2);
  // TODO: Change to Matrix3 type
  std::vector<std::vector<std::vector<Vector3d>>> _gradPhi,
      _velocity; // level set gradient and velocity on the corners
  std::vector<std::vector<std::vector<double>>>
      _phi; // Level set stored on the grid corners and its gradient values

  bool _isPressureSecondOrder;
};
} // namespace Ramuh

#endif