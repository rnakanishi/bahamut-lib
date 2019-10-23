#include <omp.h>
#include <geometry/dual_marching.h>
#include <surface/dual_cubes.h>
#include <algorithm>

namespace Leviathan {

DualCubes::DualCubes(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain)
    : LevelSetFluid3(gridSize, domain) {
  _faceSurfacePositionId = newFaceArrayLabel("facePosition");
  _faceSurfaceNormalId = newFaceArrayLabel("surfaceNormal");

  newFaceArrayLabel("facePosition");
  newFaceArrayLabel("faceNormal");
  surface = Ramuh::DualMarching3(_gridSize);
}

DualCubes::DualCubes(Leviathan::LevelSetFluid3 levelset)
    : DualCubes(levelset.getGridSize(), levelset.getDomain()) {
  swapLevelSet(levelset);
}

void DualCubes::swapLevelSet(LevelSetFluid3 levelset) {
  auto &phi = getCellScalarData("phi");
  auto &gradient = getCellArrayData("cellGradient");
  auto &lphi = levelset.getCellScalarData("phi");
  auto &lgradient = levelset.getCellArrayData("cellGradient");

  phi.clear();
  gradient.clear();
  std::copy(lphi.begin(), lphi.end(), std::back_inserter(phi));
  std::copy(lgradient.begin(), lgradient.end(), std::back_inserter(gradient));
}

void DualCubes::computeIntersectionAndNormals() {
  auto &phi = getCellScalarData("phi");
  auto &gradient = getCellArrayData("cellGradient");
  auto &ufaceLocation = getFaceArrayData(0, "facePosition");
  auto &vfaceLocation = getFaceArrayData(1, "facePosition");
  auto &wfaceLocation = getFaceArrayData(2, "facePosition");
  auto &ufaceNormals = getFaceArrayData(0, "faceNormal");
  auto &vfaceNormals = getFaceArrayData(1, "faceNormal");
  auto &wfaceNormals = getFaceArrayData(2, "faceNormal");
  auto _h = getH();

  findSurfaceCells(_h[0] * 3);
  Eigen::Array3d domainMin = _domain.getMin();
  // #pragma omp parallel for
  for (size_t surfId = 0; surfId < _surfaceCellIds.size(); surfId++) {
    int cellId = _surfaceCellIds[surfId];
    auto ijk = idToijk(cellId);
    int i = ijk[0], j = ijk[1], k = ijk[2];
    int centerSign = (phi[ijkToid(i, j, k)] < 0) ? -1 : 1;
    int id = ijkToid(i, j, k);
    if (i < _gridSize[0] - 1) {
      int neighSign = (phi[ijkToid(i + 1, j, k)] < 0) ? -1 : 1;
      if (centerSign != neighSign) {
        // Compute intersection location
        double theta = phi[ijkToid(i + 1, j, k)] - phi[ijkToid(i, j, k)];
        double x = domainMin[0] - _h[0] * phi[ijkToid(i, j, k)] / (theta) +
                   (i + 0.5) * _h[0];
        double y = domainMin[1] + (j + 0.5) * _h[1];
        double z = domainMin[2] + (k + 0.5) * _h[2];
        ufaceLocation[faceijkToid(0, i + 1, j, k)] = Eigen::Array3d(x, y, z);
        ufaceNormals[faceijkToid(0, i + 1, j, k)] =
            (gradient[id] - gradient[ijkToid(i + 1, j, k)]) * phi[id] / theta +
            gradient[id];
      }
    }
    if (j < _gridSize[1] - 1) {
      int neighSign = (phi[ijkToid(i, j + 1, k)] < 0) ? -1 : 1;
      if (centerSign != neighSign) {
        // Compute intersection location
        double theta = phi[ijkToid(i, j + 1, k)] - phi[ijkToid(i, j, k)];
        double x = domainMin[0] + (i + 0.5) * _h[0];
        double y = domainMin[1] - _h[1] * phi[ijkToid(i, j, k)] / (theta) +
                   (j + 0.5) * _h[1];
        double z = domainMin[2] + (k + 0.5) * _h[2];
        vfaceLocation[faceijkToid(1, i, j + 1, k)] = Eigen::Array3d(x, y, z);
        vfaceNormals[faceijkToid(1, i, j + 1, k)] =
            (gradient[id] - gradient[ijkToid(i, j + 1, k)]) * phi[id] / theta +
            gradient[id];
      }
    }
    if (k < _gridSize[2] - 1) {
      int neighSign = (phi[ijkToid(i, j, k + 1)] < 0) ? -1 : 1;
      if (centerSign != neighSign) {
        // Compute intersection location
        double theta = phi[ijkToid(i, j, k + 1)] - phi[ijkToid(i, j, k)];
        double x = domainMin[0] + (i + 0.5) * _h[0];
        double y = domainMin[1] + (j + 0.5) * _h[1];
        double z = domainMin[2] + -_h[2] * phi[ijkToid(i, j, k)] / (theta) +
                   (k + 0.5) * _h[2];
        wfaceLocation[faceijkToid(2, i, j, k + 1)] = Eigen::Array3d(x, y, z);
        wfaceNormals[faceijkToid(2, i, j, k + 1)] =
            (gradient[id] - gradient[ijkToid(i, j, k + 1)]) * phi[id] / theta +
            gradient[id];
      }
    }
  }
}

void DualCubes::extractSurface() {
  auto &_phi = getCellScalarData("phi");
  auto &cellGradient = getCellArrayData("cellGradient");
  auto &_ufaceLocation = getFaceArrayData(0, "facePosition");
  auto &_vfaceLocation = getFaceArrayData(1, "facePosition");
  auto &_wfaceLocation = getFaceArrayData(2, "facePosition");
  auto &_ufaceNormals = getFaceArrayData(0, "faceNormal");
  auto &_vfaceNormals = getFaceArrayData(1, "faceNormal");
  auto &_wfaceNormals = getFaceArrayData(2, "faceNormal");
  auto _h = getH();

  findSurfaceCells(_h[0] * 3);

  std::vector<std::pair<int, int>> connections;
#pragma omp parallel
  {
    Ramuh::DualMarching3 threadSurface(_gridSize);
    std::vector<std::pair<int, int>> threadConnections;
#pragma omp for
    for (size_t surfId = 0; surfId < _surfaceCellIds.size(); surfId++) {
      int cellId = _surfaceCellIds[surfId];
      auto ijk = idToijk(cellId);
      int i = ijk[0], j = ijk[1], k = ijk[2];

      std::vector<Eigen::Array3d> normalPosition;
      std::vector<Eigen::Vector3d> normal;
      bool isSurface = false;
      int id = ijkToid(i, j, k);

      // yz-plane
      if (hasSignChange(_phi[ijkToid(i, j, k)], _phi[ijkToid(i + 1, j, k)])) {
        isSurface = true;
        normal.emplace_back(_ufaceNormals[faceijkToid(0, i + 1, j, k)]);
        normalPosition.emplace_back(
            _ufaceLocation[faceijkToid(0, i + 1, j, k)]);

        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j, k), ijkToid(i, j, k - 1)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j, k), ijkToid(i, j - 1, k)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j, k - 1), ijkToid(i, j - 1, k - 1)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j - 1, k), ijkToid(i, j - 1, k - 1)));
      }
      if (hasSignChange(_phi[ijkToid(i + 1, j + 1, k)],
                        _phi[ijkToid(i, j + 1, k)])) {
        isSurface = true;
        normal.emplace_back(_ufaceNormals[faceijkToid(0, i + 1, j + 1, k)]);
        normalPosition.emplace_back(
            _ufaceLocation[faceijkToid(0, i + 1, j + 1, k)]);
      }
      if (hasSignChange(_phi[ijkToid(i, j, k + 1)],
                        _phi[ijkToid(i + 1, j, k + 1)])) {
        isSurface = true;
        normal.emplace_back(_ufaceNormals[faceijkToid(0, i + 1, j, k + 1)]);
        normalPosition.emplace_back(
            _ufaceLocation[faceijkToid(0, i + 1, j, k + 1)]);
      }
      if (hasSignChange(_phi[ijkToid(i + 1, j + 1, k + 1)],
                        _phi[ijkToid(i, j + 1, k + 1)])) {
        isSurface = true;
        normal.emplace_back(_ufaceNormals[faceijkToid(0, i + 1, j + 1, k + 1)]);
        normalPosition.emplace_back(
            _ufaceLocation[faceijkToid(0, i + 1, j + 1, k + 1)]);
      }

      // xz-plane
      if (hasSignChange(_phi[ijkToid(i + 1, j, k)],
                        _phi[ijkToid(i + 1, j + 1, k)])) {
        isSurface = true;
        normal.emplace_back(_vfaceNormals[faceijkToid(1, i + 1, j + 1, k)]);
        normalPosition.emplace_back(
            _vfaceLocation[faceijkToid(1, i + 1, j + 1, k)]);
      }
      if (hasSignChange(_phi[ijkToid(i, j + 1, k)], _phi[ijkToid(i, j, k)])) {
        isSurface = true;
        normal.emplace_back(_vfaceNormals[faceijkToid(1, i, j + 1, k)]);
        normalPosition.emplace_back(
            _vfaceLocation[faceijkToid(1, i, j + 1, k)]);

        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j, k), ijkToid(i, j, k - 1)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j, k), ijkToid(i - 1, j, k)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j, k - 1), ijkToid(i - 1, j, k - 1)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i - 1, j, k), ijkToid(i - 1, j, k - 1)));
      }
      if (hasSignChange(_phi[ijkToid(i + 1, j, k + 1)],
                        _phi[ijkToid(i + 1, j + 1, k + 1)])) {
        isSurface = true;
        normal.emplace_back(_vfaceNormals[faceijkToid(1, i + 1, j + 1, k + 1)]);
        normalPosition.emplace_back(
            _vfaceLocation[faceijkToid(1, i + 1, j + 1, k + 1)]);
      }
      if (hasSignChange(_phi[ijkToid(i, j + 1, k + 1)],
                        _phi[ijkToid(i, j, k + 1)])) {
        isSurface = true;
        normal.emplace_back(_vfaceNormals[faceijkToid(1, i, j + 1, k + 1)]);
        normalPosition.emplace_back(
            _vfaceLocation[faceijkToid(1, i, j + 1, k + 1)]);
      }

      // xy-plane
      if (hasSignChange(_phi[ijkToid(i, j, k)], _phi[ijkToid(i, j, k + 1)])) {
        isSurface = true;
        normal.emplace_back(_wfaceNormals[faceijkToid(2, i, j, k + 1)]);
        normalPosition.emplace_back(
            _wfaceLocation[faceijkToid(2, i, j, k + 1)]);

        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j, k), ijkToid(i - 1, j, k)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j, k), ijkToid(i, j - 1, k)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i - 1, j, k), ijkToid(i - 1, j - 1, k)));
        threadConnections.emplace_back(
            std::make_pair(ijkToid(i, j - 1, k), ijkToid(i - 1, j - 1, k)));
      }
      if (hasSignChange(_phi[ijkToid(i + 1, j, k)],
                        _phi[ijkToid(i + 1, j, k + 1)])) {
        isSurface = true;
        normal.emplace_back(_wfaceNormals[faceijkToid(2, i + 1, j, k + 1)]);
        normalPosition.emplace_back(
            _wfaceLocation[faceijkToid(2, i + 1, j, k + 1)]);
      }
      if (hasSignChange(_phi[ijkToid(i + 1, j + 1, k)],
                        _phi[ijkToid(i + 1, j + 1, k + 1)])) {
        isSurface = true;
        normal.emplace_back(_wfaceNormals[faceijkToid(2, i + 1, j + 1, k + 1)]);
        normalPosition.emplace_back(
            _wfaceLocation[faceijkToid(2, i + 1, j + 1, k + 1)]);
      }
      if (hasSignChange(_phi[ijkToid(i, j + 1, k)],
                        _phi[ijkToid(i, j + 1, k + 1)])) {
        isSurface = true;
        normal.emplace_back(_wfaceNormals[faceijkToid(2, i, j + 1, k + 1)]);
        normalPosition.emplace_back(
            _wfaceLocation[faceijkToid(2, i, j + 1, k + 1)]);
      }
      // Solve quadratic function
      if (isSurface) {
        Eigen::Array3d cubeMin =
            _domain.getMin() +
            _h.cwiseProduct(Eigen::Array3d(i + 0.5, j + 0.5, k + 0.5));
        Eigen::Array3d cubeMax =
            _domain.getMin() +
            _h.cwiseProduct(Eigen::Array3d(i + 1.5, j + 1.5, k + 1.5));
        int thread = omp_get_thread_num();
        auto x = threadSurface.evaluateCube(
            Eigen::Array3i(i, j, k), normalPosition, normal,
            Ramuh::BoundingBox3(cubeMin, cubeMax));
        // std::cout << x[0] << " " << x[1] << " " << x[2] <<
        // std::endl;
      }
    }
#pragma omp critical
    {
      surface.merge(threadSurface);
      connections.insert(connections.end(), threadConnections.begin(),
                         threadConnections.end());
    }
  }
  surface.setBaseFolder(_baseFolder);
  surface.reconstruct(connections);
}

void DualCubes::setFolder(std ::string folder) { _baseFolder = folder; }

void DualCubes::resetFileCounter() { surface.resetCounter(); }

bool DualCubes::hasSignChange(double valueA, double valueB) {
  int signA = (valueA < 0) ? -1 : 1;
  int signB = (valueB < 0) ? -1 : 1;
  return signA != signB;
}

} // namespace Leviathan