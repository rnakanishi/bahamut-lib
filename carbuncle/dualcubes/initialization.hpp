
#include <surface/dual_squares.h>
#include <structures/particle_system2.h>

enum class ParametricSurface : int { CIRCLE, SQUARE };

void defineCellsVelocity(Leviathan::DualSquares &levelset) {
  auto &u = levelset.getFaceScalarData(0, "faceVelocity");
  for (size_t faceId = 0; faceId < levelset.faceCount(0); faceId++) {
    auto ij = levelset.faceIdToij(0, faceId);
    auto p = levelset.getFacePosition(0, faceId);
    u[faceId] = -p[1];
  }
  auto &v = levelset.getFaceScalarData(1, "faceVelocity");
  for (size_t faceId = 0; faceId < levelset.faceCount(1); faceId++) {
    auto ij = levelset.faceIdToij(1, faceId);
    auto p = levelset.getFacePosition(1, faceId);
    v[faceId] = p[0];
  }
}

void defineParticlesVelocity(Ramuh::ParticleSystem2 &particles) {
  auto &velocity = particles.getParticleArrayData("particleVelocity");
  for (size_t pid = 0; pid < particles.particleCount(); pid++) {
    auto p = particles.getParticlePosition(pid);
    velocity[pid] = Eigen::Array2d(-p[1], p[0]);
  }
}

void initializeCube(Leviathan::DualSquares &squares, Eigen::Array2d center,
                    double radius, ParametricSurface surface) {

  Eigen::Array2d domainMin = squares.getDomain().getMin();
  auto &_phi = squares.getCellScalarData("phi");
  auto gridSize = squares.getGridSize();
  auto h = squares.getH();

  // ====== Initialize cell surfaces
  for (int j = 0; j < gridSize[1]; j++)
    for (int i = 0; i < gridSize[0]; i++) {
      Eigen::Array2d position =
          domainMin + Eigen::Array2d(i, j).cwiseProduct(h) + h / 2.0;
      // double distance = (position - center).matrix().norm() - radius;

      double x, y, x2, y2;
      double distance;
      x = position[0] - center[0];
      y = position[1] - center[1];
      x2 = x * x;
      y2 = y * y;

      double r2;
      double r;
      switch (surface) {
      // SPHERE
      case ParametricSurface::CIRCLE:
        distance = x2 + y2 - radius * radius;
        break;
      case ParametricSurface::SQUARE:
        // CUBE
        distance = std::max(std::fabs(x), std::fabs(y)) - radius;
        if (distance > 0) {
          position = position.abs();
          distance = 0.0;
          x = std::max(0.0, position[0] - radius);
          y = std::max(0.0, position[1] - radius);
          distance = sqrt(x * x + y * y);
        }
        break;
      default:
        distance = 1e8;
      }
      _phi[squares.ijToid(i, j)] =
          std::min(_phi[squares.ijToid(i, j)], distance);
    }
}

void initializeGradientsAtIntersection(Leviathan::DualSquares &squares,
                                       Eigen::Array2d center, double radius,
                                       ParametricSurface surface) {
  auto &phi = squares.getCellScalarData("phi");
  auto &gradient = squares.getCellArrayData("cellGradient");
  auto &ufaceLocation = squares.getFaceArrayData(0, "facePosition");
  auto &vfaceLocation = squares.getFaceArrayData(1, "facePosition");
  auto &ufaceNormals = squares.getFaceArrayData(0, "faceNormal");
  auto &vfaceNormals = squares.getFaceArrayData(1, "faceNormal");
  auto _h = squares.getH();
  auto gridSize = squares.getGridSize();

  auto surfaceCellIds = squares.findSurfaceCells(_h[0] * 8);
  Eigen::Array2d domainMin = squares.getDomain().getMin();
  // #pragma omp parallel for
  // for (size_t cellId = 0; cellId < cellCount(); cellId++) {
  for (size_t surfId = 0; surfId < surfaceCellIds.size(); surfId++) {
    int cellId = surfaceCellIds[surfId];
    auto ij = squares.idToij(cellId);
    int i = ij[0], j = ij[1];
    int centerSign = (phi[squares.ijToid(i, j)] < 0) ? -1 : 1;
    int id = squares.ijToid(i, j);
    if (i < gridSize[0] - 1) {
      if (squares.hasSignChange(phi[cellId], phi[squares.ijToid(i + 1, j)])) {
        // Compute intersection location
        double theta =
            phi[squares.ijToid(i + 1, j)] - phi[squares.ijToid(i, j)];
        double x = domainMin[0] - _h[0] * phi[squares.ijToid(i, j)] / (theta) +
                   (i + 0.5) * _h[0];
        double y = domainMin[1] + (j + 0.5) * _h[1];
        ufaceLocation[squares.faceijToid(0, i + 1, j)] = Eigen::Array2d(x, y);

        Eigen::Array2d normal = Eigen::Array2d(0);
        if (x < center[0])
          normal[0] = -1;
        else
          normal[0] = 1;
        ufaceNormals[squares.faceijToid(0, i + 1, j)] = normal;
        int ii = 0;
      }
    }

    if (j < gridSize[1] - 1) {
      if (squares.hasSignChange(phi[cellId], phi[squares.ijToid(i, j + 1)])) {
        // Compute intersection location
        double theta =
            phi[squares.ijToid(i, j + 1)] - phi[squares.ijToid(i, j)];
        double x = domainMin[0] + (i + 0.5) * _h[0];
        double y = domainMin[1] - _h[1] * phi[squares.ijToid(i, j)] / (theta) +
                   (j + 0.5) * _h[1];
        vfaceLocation[squares.faceijToid(1, i, j + 1)] = Eigen::Array2d(x, y);

        Eigen::Array2d normal = Eigen::Array2d(0);
        if (y < center[1])
          normal[1] = -1;
        else
          normal[1] = 1;
        vfaceNormals[squares.faceijToid(1, i, j + 1)] = normal;
        int ii = 0;
      }
    }

    // Computing cell center gradient
    Eigen::Array2d normalDirection(0);
    auto p = squares.getCellPosition(cellId);
    p -= center;

    int signal = 1;
    if (std::abs(p[0]) < radius && std::abs(p[1]) < radius) {
      // Inside cube
      if (std::abs(p[0]) > std::abs(p[1])) {
        signal = (p[0]) < 0 ? -1 : 1;
        normalDirection = Eigen::Array2d(1, 0);
      } else {
        signal = (p[1]) < 0 ? -1 : 1;
        normalDirection = Eigen::Array2d(0, 1);
      }
    } else {
      // Diagonal values take direction of the nearest face
      if (std::abs(p[0]) < std::abs(p[1])) {
        signal = (p[1]) < 0 ? -1 : 1;
        normalDirection = Eigen::Array2d(0, 1);
      } else {
        signal = (p[0]) < 0 ? -1 : 1;
        normalDirection = Eigen::Array2d(1, 0);
      }
    }
    gradient[cellId] = normalDirection * signal;
  }
}

void createLineMesh(Ramuh::LineMesh &mesh, Eigen::Array2d center,
                    double radius) {

  // Creating a single square outline
  Ramuh::BoundingBox2 bbox(center - radius, center + radius);
  Eigen::Array2d min, max;
  min = bbox.getMin();
  max = bbox.getMax();
  std::vector<int> points(4);

  points[0] = mesh.addVertex(min);
  points[1] = mesh.addVertex(Eigen::Array2d(max[0], min[1]));
  points[2] = mesh.addVertex(max);
  points[3] = mesh.addVertex(Eigen::Array2d(min[0], max[1]));

  mesh.connectVertices(Eigen::Array2i(points[0], points[1]));
  mesh.connectVertices(Eigen::Array2i(points[1], points[2]));
  mesh.connectVertices(Eigen::Array2i(points[2], points[3]));
  mesh.connectVertices(Eigen::Array2i(points[3], points[0]));
}
