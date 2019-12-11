#include <omp.h>
#include <geometry/dual_marching.h>
#include <surface/dual_squares.h>
#include <algorithm>

namespace Leviathan {

DualSquares::DualSquares(Eigen::Array2i gridSize, Ramuh::BoundingBox2 domain)
    : LevelSetFluid2(gridSize, domain) {
  _faceSurfacePositionId = newFaceArrayLabel("facePosition");
  _faceSurfaceNormalId = newFaceArrayLabel("surfaceNormal");

  newFaceArrayLabel("facePosition");
  newFaceArrayLabel("faceNormal");
  surface = Ramuh::DualMarching2(_gridSize);
}

DualSquares::DualSquares(Leviathan::LevelSetFluid2 levelset)
    : DualSquares(levelset.getGridSize(), levelset.getDomain()) {
  swapLevelSet(levelset);
}

void DualSquares::swapLevelSet(LevelSetFluid2 &levelset) {
  auto &phi = getCellScalarData("phi");
  auto &gradient = getCellArrayData("cellGradient");
  auto &lphi = levelset.getCellScalarData("phi");
  auto &lgradient = levelset.getCellArrayData("cellGradient");

  phi.clear();
  gradient.clear();
  std::copy(lphi.begin(), lphi.end(), std::back_inserter(phi));
  std::copy(lgradient.begin(), lgradient.end(), std::back_inserter(gradient));

  surface = Ramuh::DualMarching2(_gridSize);
}

void DualSquares::computeIntersectionAndNormals() {
  auto &phi = getCellScalarData("phi");
  auto &gradient = getCellArrayData("cellGradient");
  auto &ufaceLocation = getFaceArrayData(0, "facePosition");
  auto &vfaceLocation = getFaceArrayData(1, "facePosition");
  auto &ufaceNormals = getFaceArrayData(0, "faceNormal");
  auto &vfaceNormals = getFaceArrayData(1, "faceNormal");
  auto _h = getH();

  findSurfaceCells(_h[0] * 8);
  Eigen::Array2d domainMin = _domain.getMin();
  // #pragma omp parallel for
  // for (size_t cellId = 0; cellId < cellCount(); cellId++) {
  for (size_t surfId = 0; surfId < _surfaceCellIds.size(); surfId++) {
    int cellId = _surfaceCellIds[surfId];
    auto ij = idToij(cellId);
    int i = ij[0], j = ij[1];
    int centerSign = (phi[ijToid(i, j)] < 0) ? -1 : 1;
    int id = ijToid(i, j);
    if (i < _gridSize[0] - 1) {
      if (hasSignChange(phi[cellId], phi[ijToid(i + 1, j)])) {
        // Compute intersection location
        double theta = phi[ijToid(i + 1, j)] - phi[ijToid(i, j)];
        double x = domainMin[0] - _h[0] * phi[ijToid(i, j)] / (theta) +
                   (i + 0.5) * _h[0];
        double y = domainMin[1] + (j + 0.5) * _h[1];
        ufaceLocation[faceijToid(0, i + 1, j)] = Eigen::Array2d(x, y);
        ufaceNormals[faceijToid(0, i + 1, j)] =
            (gradient[id] - gradient[ijToid(i + 1, j)]) * phi[id] / theta +
            gradient[id];
        int ii = 0;
      }
    }
    if (j < _gridSize[1] - 1) {
      if (hasSignChange(phi[cellId], phi[ijToid(i, j + 1)])) {
        // Compute intersection location
        double theta = phi[ijToid(i, j + 1)] - phi[ijToid(i, j)];
        double x = domainMin[0] + (i + 0.5) * _h[0];
        double y = domainMin[1] - _h[1] * phi[ijToid(i, j)] / (theta) +
                   (j + 0.5) * _h[1];
        vfaceLocation[faceijToid(1, i, j + 1)] = Eigen::Array2d(x, y);
        vfaceNormals[faceijToid(1, i, j + 1)] =
            (gradient[id] - gradient[ijToid(i, j + 1)]) * phi[id] / theta +
            gradient[id];
        int ii = 0;
      }
    }
  }
}

void DualSquares::extractSurface() {
  auto &_phi = getCellScalarData("phi");
  auto &cellGradient = getCellArrayData("cellGradient");
  auto &_ufaceLocation = getFaceArrayData(0, "facePosition");
  auto &_vfaceLocation = getFaceArrayData(1, "facePosition");
  auto &_ufaceNormals = getFaceArrayData(0, "faceNormal");
  auto &_vfaceNormals = getFaceArrayData(1, "faceNormal");
  auto _h = getH();

  surface.clear();
  findSurfaceCells(_h[0] * 8);

  std::vector<std::pair<int, int>> connections;
#pragma omp parallel
  {
    Ramuh::DualMarching2 threadSurface(_gridSize);
    std::vector<std::pair<int, int>> threadConnections;
#pragma omp for
    for (size_t cellId = 0; cellId < cellCount(); cellId++) {
      // The resulting point is stored in the evaluating cell
      // Connections between cells are made with those that share an
      // intersecting edge
      auto ij = idToij(cellId);
      int i = ij[0], j = ij[1];

      std::vector<Eigen::Array2d> normalPosition;
      std::vector<Eigen::Vector2d> normal;
      bool isSurface = false;
      int id = ijToid(i, j);
      if (!(i + 1 < _gridSize[0] && j + 1 < _gridSize[1]))
        continue;

      if ((i + 1 < _gridSize[0]) &&
          hasSignChange(_phi[ijToid(i, j)], _phi[ijToid(i + 1, j)])) {
        isSurface = true;
        normal.emplace_back(_ufaceNormals[faceijToid(0, i + 1, j)]);
        normalPosition.emplace_back(_ufaceLocation[faceijToid(0, i + 1, j)]);
      }
      if ((i + 1 < _gridSize[0] && j + 1 < _gridSize[1]) &&
          hasSignChange(_phi[ijToid(i + 1, j + 1)], _phi[ijToid(i, j + 1)])) {
        isSurface = true;
        normal.emplace_back(_ufaceNormals[faceijToid(0, i + 1, j + 1)]);
        normalPosition.emplace_back(
            _ufaceLocation[faceijToid(0, i + 1, j + 1)]);
        // Connections are made between cells that share an edge containing
        // intersection
        threadConnections.emplace_back(
            std::make_pair(ijToid(i, j + 1), ijToid(i, j)));
      }

      if ((i + 1 < _gridSize[0] && j + 1 < _gridSize[1]) &&
          hasSignChange(_phi[ijToid(i + 1, j + 1)], _phi[ijToid(i + 1, j)])) {
        isSurface = true;
        normal.emplace_back(_vfaceNormals[faceijToid(1, i + 1, j + 1)]);
        normalPosition.emplace_back(
            _vfaceLocation[faceijToid(1, i + 1, j + 1)]);
        threadConnections.emplace_back(
            std::make_pair(ijToid(i, j), ijToid(i + 1, j)));
      }
      if ((j + 1 < _gridSize[1]) &&
          hasSignChange(_phi[ijToid(i, j)], _phi[ijToid(i, j + 1)])) {
        isSurface = true;
        normal.emplace_back(_vfaceNormals[faceijToid(1, i, j + 1)]);
        normalPosition.emplace_back(_vfaceLocation[faceijToid(1, i, j + 1)]);
      }

      // Solve quadratic function
      if (isSurface) {
        Eigen::Array2d cubeMin =
            _domain.getMin() +
            _h.cwiseProduct(Eigen::Array2d(i + 0.5, j + 0.5));
        Eigen::Array2d cubeMax =
            _domain.getMin() +
            _h.cwiseProduct(Eigen::Array2d(i + 1.5, j + 1.5));
        int thread = omp_get_thread_num();
        auto x = threadSurface.evaluateSquare(
            Eigen::Array2i(i, j), normalPosition, normal,
            Ramuh::BoundingBox2(cubeMin, cubeMax));
      }
    }
#pragma omp critical(surfaceMerge)
    {
      surface.merge(threadSurface);
      connections.insert(connections.end(), threadConnections.begin(),
                         threadConnections.end());
    }
  }
  surface.setBaseFolder(_baseFolder);
  surface.reconstruct(connections);
}

void DualSquares::setFolder(std::string folder) { _baseFolder = folder; }

void DualSquares::resetFileCounter() { surface.resetCounter(); }

void DualSquares::resetFileCounter(int value) { surface.resetCounter(value); }

bool DualSquares::hasSignChange(double valueA, double valueB) {
  int signA = (valueA < 0) ? -1 : 1;
  int signB = (valueB < 0) ? -1 : 1;
  return signA != signB;
}

} // namespace Leviathan