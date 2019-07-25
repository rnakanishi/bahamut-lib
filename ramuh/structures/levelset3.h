#ifndef __RAMUH_LEVELSET3_H__
#define __RAMUH_LEVELSET3_H__

#include <structures/grid3.h>
#include <structures/triangle_mesh.h>
#include <utils/material.h>
#include <Eigen/Dense>
#include <glm/glm.hpp>

namespace Ramuh {

class LevelSet3 : public RegularGrid3 {
public:
  LevelSet3();
  LevelSet3(Eigen::Array3i resolution);

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
  void addSphereSurface(Eigen::Array3d center, double radius);
  void addCubeSurface(Eigen::Array3d lower, Eigen::Array3d upper);

  /**
   * @brief Perform marching tetrahedra algorithm looking for zero level set.
   * The coordinates of the vertices correspond to the cell center. Possibly,
   * created vertices may be doubled and faces completely independent
   *
   * @return TriangleMesh full model with vertices and faces of the zero
   *levelset
   **/
  TriangleMesh marchingTetrahedra();

  /**
   * @brief Set the Pressure Second Order value. If true, then second order for
   *free surface will be used
   *
   * @param value
   **/
  void setPressureSecondOrder(bool value);

  /**
   * @brief For each cell, evaluates levelset values. If a levelset is zero or
   *negative, then it is set as fluid cell.
   *
   **/
  void checkCellMaterial();

  void printVertexVelocity();

  void writeLevelSetValue(const char *filename);
  void readLevelSetValue(const char *filename);

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
  void _triangulate(glm::vec3 coords[4], double values[4], TriangleMesh &mesh);

  /**
   * @brief Interpolate phi values using proper cubic interpolator. This method
   *builds the stencil for interpolation based on @position value.
   *
   * @param position Location where the value for phi is wanted
   * @return double interpolated value
   **/
  double _interpolatePhi(Eigen::Array3d position);
  double _interpolatePhi(Eigen::Array3d position, double &min, double &max);
  double _interpolatePhi(Eigen::Array3d position, int signal);

  glm::vec3 _findSurfaceCoordinate(glm::ivec3 v1, glm::ivec3 v2);
  Eigen::Array3d _findSurfaceCoordinate(Eigen::Array3i v1, Eigen::Array3i v2);
  glm::vec3 _findSurfaceCoordinate(glm::vec3 coord1, glm::vec3 coord2,
                                   float phi1, float phi2);

  // TODO: Change to Matrix3 type
  std::vector<std::vector<std::vector<Vector3d>>>
      _velocity; // level set gradient and velocity on the corners
  std::vector<std::vector<std::vector<double>>>
      _phi; // Level set stored on the grid corners and its gradient values

  bool _isPressureSecondOrder;
};
} // namespace Ramuh

#endif