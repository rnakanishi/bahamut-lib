#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <iomanip>
#include "dual_cubes.h"
#include <surface/dual_squares.h>
#include <structures/line_mesh.hpp>
#include <fstream>
#include "initialization.hpp"
#include "tangent_particles2.hpp"
#include "levelset_from_mesh.hpp"
#include <Eigen/Geometry>

void printParticles(Carbuncle::TangentParticles2 particles) {
  std::ofstream file("matlab/particles.txt");
  auto &normalPair = particles.getPairMap();
  for (auto normal : normalPair) {
    if (particles.isActive(normal.first)) {
      auto origin = particles.getParticlePosition(normal.first);
      auto ending = particles.getParticlePosition(normal.second);
      Eigen::Array2d tangent((ending - origin).matrix());
      Eigen::Vector2d direction(tangent[1], -tangent[0]);
      file << origin.transpose() << " " << direction.normalized().transpose()
           << std::endl;
    }
  }
  file.close();
}

void printParticlesPerCell(Carbuncle::TangentParticles2 particles) {
  std::ofstream file("matlab/particles.txt");
  auto &normalPair = particles.getPairMap();
  auto cells = particles.computeCellsWithParticles();
  for (auto cell : cells) {
    auto pVector = particles.getParticlesInCell(cell);
    for (auto particleId : pVector) {
      if (normalPair.find(particleId) == normalPair.end())
        continue;
      auto normal = normalPair[particleId];
      if (particles.isActive(particleId) && particles.isActive(normal)) {
        auto origin = particles.getParticlePosition(particleId);
        auto ending = particles.getParticlePosition(normal);
        Eigen::Array2d tangent((ending - origin).matrix());
        Eigen::Vector2d direction(tangent[1], -tangent[0]);
        file << origin.transpose() << " " << direction.normalized().transpose()
             << std::endl;
      }
    }
  }
  file.close();
}

void fixLevelsetFromMesh(Leviathan::DualSquares &levelset,
                         Ramuh::LineMesh &mesh) {
  auto &phi = levelset.getCellScalarData("phi");

  // For all segments
  auto &segments = mesh.getSegmentsList();
  auto &vertices = mesh.getVerticesList();
  Eigen::Vector2d tangent;
  auto h = levelset.getH();

  for (auto segment : segments) {
    // Get the origin point positions and find which cell it belongs
    Eigen::Array2d origin = mesh.getVertexPosition(segment[0]);
    Eigen::Array2d ending = mesh.getVertexPosition(segment[1]);
    double angle = (ending[1] - origin[1]) / (ending[0] - origin[0]);
    auto cellId = levelset.findCellIdByCoordinate(origin);
    auto cellij = levelset.idToij(cellId);

    // Do rasterization for each cell the segment crosses
    bool ended = false;
    Eigen::Array2d newOrigin = origin, newEnding;
    Ramuh::BoundingBox2 cellBbox = levelset.getCellBoundingBox(cellId);
    std::vector<bool> visited(phi.size(), false);
    while (!ended) {
      Eigen::Vector2i cellInc(0, 0);
      newEnding = ending;
      // find which square edge intersects the surface
      if (!cellBbox.contains(ending)) {
        auto cellMin = cellBbox.getMin();
        auto cellMax = cellBbox.getMax();
        // for each edge, check intersection. If true, project to that edge:
        // bottom/top
        if (newEnding[1] < cellMin[1]) {
          newEnding[1] = cellMin[1];
          newEnding[0] = (newEnding[1] - newOrigin[1]) / angle + origin[0];
          if (cellBbox.contains(newEnding))
            cellij[1]--;
        } else if (newEnding[1] > cellMax[1]) {
          newEnding[1] = cellMax[1];
          newEnding[0] = (newEnding[1] - newOrigin[1]) / angle + origin[0];
          if (cellBbox.contains(newEnding))
            cellij[1]++;
        }
        // left/right
        if (newEnding[0] < cellMin[0]) {
          newEnding[0] = cellMin[0];
          newEnding[1] = angle * (newEnding[0] - newOrigin[0]) + newOrigin[1];
          if (cellBbox.contains(newEnding))
            cellij[0]--;
        } else if (newEnding[0] > cellMax[0]) {
          newEnding[0] = cellMax[0];
          newEnding[1] = angle * (newEnding[0] - newOrigin[0]) + newOrigin[1];
          if (cellBbox.contains(newEnding))
            cellij[0]++;
        }
      } else {
        ended = true;
      }

      // Compute cell center distance to the segment
      auto neighborCells = levelset.findCellNeighbors(cellId, 2);
      for (auto cell : neighborCells) {
        Eigen::Array2d cellCenter = levelset.getCellPosition(cell);
        Eigen::Vector2d direction, target;
        direction = newEnding - newOrigin;
        target = cellCenter - newOrigin;
        int pSignal = (phi[cell] < 0) ? -1 : 1;
        double distance;
        double cross = (target[0] * direction[1] - target[1] * direction[0]);
        int dSignal;
        if (cross == 0) {
          distance = 0.0;
          visited[cell] = true;
          phi[cell] = 0.0;
        } else {
          if (cross * pSignal < 0) {
            distance = distanceToSegment(newOrigin, newEnding, cellCenter);
          } else {
            distance = distanceToSegment(newEnding, newOrigin, cellCenter);
          }
          dSignal = (distance < 0) ? -1 : 1;
          if (dSignal == pSignal)
            if (!visited[cell]) {
              visited[cell] = true;
              phi[cell] = distance;
            } else {
              phi[cell] =
                  std::min(std::abs(phi[cell]), std::abs(distance)) * pSignal;
            }
        }
      }
      // Special case has to be treated when segment has reach an end: check for
      // singularity
      if (ended) {
      }

      // Prepare for next cell or segment
      cellId = levelset.ijToid(cellij[0], cellij[1]);
      cellBbox = levelset.getCellBoundingBox(cellId);
      newOrigin = newEnding;
    }
  }
}

void writeMesh(Ramuh::LineMesh &mesh, std::string baseFolder, int count) {

  std::ofstream file;
  std::stringstream filename;
  filename << baseFolder << std::setfill('0') << std::setw(4) << count++
           << ".obj";
  auto fullpath = filename.str();
  file.open(fullpath.c_str(), std::ofstream::out);
  if (!file.is_open()) {
    std::cerr << "\033[1;31m[ X ]\033[0m Failed opening file " << filename.str()
              << std::endl;
    return;
  }

  auto &vertices = mesh.getVerticesList();
  auto &segments = mesh.getSegmentsList();
  int nVertices = vertices.size();
  std::vector<int> inactiveVertices;
  inactiveVertices.clear();
  for (size_t vId = 0; vId < vertices.size(); vId++) {
    auto vertex = vertices[vId];
    if (mesh.isVertexActive(vId))
      file << "v " << vertex[0] << " " << vertex[1] << " 0" << std::endl;
    else {
      inactiveVertices.emplace_back(vId);
      nVertices--;
    }
  }
  for (size_t vId = 0; vId < vertices.size(); vId++) {
    auto vertex = vertices[vId];
    if (mesh.isVertexActive(vId))
      file << "v " << vertex[0] << " " << vertex[1] << " 1" << std::endl;
  }

  for (size_t sId = 0; sId < segments.size(); sId++) {
    if (!mesh.isSegmentActive(sId))
      continue;
    auto segment = segments[sId];
    Eigen::Array2i reduction(0, 0);
    for (auto vertex : inactiveVertices) {
      if (segment[0] >= vertex)
        reduction[0]++;
      if (segment[1] >= vertex)
        reduction[1]++;
    }
    segment = segment - reduction;

    file << "f " << segment[0] + 1 << " " << segment[1] + 1 << " "
         << segment[1] + 1 + nVertices << " " << segment[0] + 1 + nVertices
         << " " << std::endl;
  }
  std::cerr << "File written: " << filename.str();
  std::cerr << ": " << vertices.size() << " vertices, " << segments.size()
            << " faces.\n";
}

void printLevelset(Leviathan::DualSquares squares) {
  std::ofstream file("matlab/levelset.txt");
  auto &phi = squares.getCellScalarData("phi");
  for (size_t cellId = 0; cellId < squares.cellCount(); cellId++) {
    file << phi[cellId] << " ";
  }

  file.close();
}

std::map<int, std::vector<int>>
findCellSegments(Leviathan::DualSquares &levelset, Ramuh::LineMesh &mesh) {
  auto segments = mesh.getSegmentsList();
  std::map<int, std::vector<int>> cellSegments;
  for (size_t segId = 0; segId < mesh.getSegmentsCount(); segId++) {
    auto vertices = segments[segId];
    auto vPos = mesh.getVertexPosition(vertices[0]);
    int cellId = levelset.findCellIdByCoordinate(vPos);
    // check if segments are added corretly
    cellSegments[cellId].emplace_back(segId);
  }
  return cellSegments;
}

void resampleParticlesFromGeometry(Carbuncle::TangentParticles2 &particles,
                                   Leviathan::DualSquares &levelset,
                                   Ramuh::LineMesh &mesh) {
  // Check particle interpolated levelset
  auto cellsWithParticles = particles.computeCellsWithParticles();
  auto phi = levelset.getCellScalarData("phi");
  auto cellSegments = findCellSegments(levelset, mesh);
  std::vector<int> cellsToAddParticles;

  for (auto cell : cellsWithParticles) {
    int cellId = cell;

    if (!levelset.isSurfaceCell(cellId)) {
      // If interior cell, remove all particles
      auto allParticles = particles.getParticlesInCell(cell);
      for (auto particleId : allParticles) {
        particles.removeParticle(particleId);
        particles.removeParticle(particles.getParticlePairId(particleId));
      }
    } else {
      // Else if surface cell, check particles over segments
      // Main problem happens when particles form a corner

      // reseed particles over the segment intersection
      auto particlesInCell = particles.getParticlesInCell(cell);
      if (particlesInCell.size() < 10) {
        // Reseed particles
        for (auto pId : particlesInCell) {
          particles.removeParticle(pId);
        }
        auto segments = cellSegments[cell];
        auto bbox = levelset.getCellBoundingBox(cell);

        for (auto segment : segments) {
          auto vertices = mesh.getSegmentVertices(segment);
          auto intersec =
              bbox.findIntersection(mesh.getVertexPosition(vertices[0]),
                                    mesh.getVertexPosition(vertices[1]));
          Eigen::Array2d inside = mesh.getVertexPosition(vertices[0]);
          if (!bbox.contains(inside))
            inside = mesh.getVertexPosition(vertices[1]);
          particles.seedParticleOverSegment(inside, intersec,
                                            levelset.getH()[0]);
        }
      }
      // Measure the segment length

      if (phi[cell] < 0) {
        auto allParticles = particles.getParticlesInCell(cell);
        auto particlePosition =
            particles.getParticleArrayData("particlePosition");
        auto h = levelset.getH();

        // Interpolate particle values
        for (auto particleId : allParticles) {
          if (particles.isActive(particleId)) {
            double particlePhi = levelset.interpolateCellScalarData(
                "phi", particlePosition[particleId]);
            if (particlePhi < -.5 * h[0])
              // std::cerr << "remove";
              particles.removeParticle(particleId);
          }
        }
      }

      // count the number of particles and compute the ratio of particles per
      // cell per segment Iff the count is below a threshold, remove all
      // partices and refill the cell using the semgmebt informations
      // Careful with corner segments
    }
  }
}

void resampleParticlesFromLevelset(Carbuncle::TangentParticles2 &particles,
                                   Leviathan::DualSquares levelset) {

  // Check particle interpolated levelset
  auto cellsWithParticles = particles.computeCellsWithParticles();
  auto phi = levelset.getCellScalarData("phi");

  for (auto cell : cellsWithParticles) {
    int cellId = cell;

    if (!levelset.isSurfaceCell(cellId)) {
      // If interior cell, remove all particles
      auto allParticles = particles.getParticlesInCell(cell);
      for (auto particleId : allParticles) {
        particles.removeParticle(particleId);
        particles.removeParticle(particles.getParticlePairId(particleId));
      }
    } else {
      // Else if surface cell, check particles over segments
      // Main problem happens when particles form a corner
      if (phi[cell] < 0) {
        auto allParticles = particles.getParticlesInCell(cell);
        auto particlePosition =
            particles.getParticleArrayData("particlePosition");
        auto h = levelset.getH();

        // Interpolate particle values
        for (auto particleId : allParticles) {
          if (particles.isActive(particleId)) {
            double particlePhi = levelset.interpolateCellScalarData(
                "phi", particlePosition[particleId]);
            if (particlePhi < -.15 * h[0])
              // std::cerr << "remove";
              particles.removeParticle(particleId);
          }
        }
      }
    }
  }
}

void removeParticlesFromSquare(Carbuncle::TangentParticles2 &particles,
                               Eigen::Array2d center, double radius,
                               double degree) {
  Eigen::Transform<double, 2, Eigen::Affine> t;
  t.setIdentity();
  t.translate(center.matrix());
  t.rotate(degree * M_PI / 180);
  t.scale(radius);
  std::cerr << t.matrix() << std::endl;

  auto &positions = particles.getParticleArrayData("particlePosition");
  auto &tanMap = particles.getPairMap();

  for (int pid = 0; pid < positions.size(); pid++) {
    auto position = positions[pid];
    Eigen::Vector2d newPos = t.inverse() * position.matrix();
    if (newPos.cwiseAbs()[0] < 1 && newPos.cwiseAbs()[1] < 1) {
      // Inner particle: remove
      particles.removeParticle(pid);
      particles.removeParticle(tanMap[pid]);
    }
  }
}

int main(int argc, char const *argv[]) {
  Leviathan::DualSquares cubes(
      Eigen::Array2i(65, 65),
      Ramuh::BoundingBox2(Eigen::Array2d(-5, -5), Eigen::Array2d(5, 5)));
  Leviathan::DualSquares cubes2(cubes);
  Carbuncle::TangentParticles2 particles(cubes);
  Carbuncle::TangentParticles2 particles2(cubes);
  std::vector<int> segmentParticles;
  Eigen::Array2d center = Eigen::Array2d(1.35, 0);
  double radius = 1.20001;

  initializeCube(cubes, center, radius, ParametricSurface::SQUARE);
  initializeCube(cubes2, -center, radius, ParametricSurface::SQUARE);
  initializeGradientsAtIntersection(cubes, center, radius,
                                    ParametricSurface::SQUARE);
  cubes.setFolder("results/dualSquares/");
  Ramuh::LineMesh squareMesh, squareMesh2, finalMesh;
  createLineMesh(squareMesh, center, radius);
  createLineMesh(squareMesh2, -center, radius);

  defineCellsVelocity(cubes);

  particles.defineRotationCenter(center);
  particles.seedParticlesOverSurface(cubes, squareMesh);
  particles2.defineRotationCenter(-center);
  particles2.seedParticlesOverSurface(cubes, squareMesh2);
  defineParticlesVelocity(particles);
  defineParticlesVelocity(particles2);
  printParticles(particles);

  // cubes.computeCellsGradient();

  squareMesh = particles.extractSurface(cubes);
  extractLevelsetFromMesh(cubes, squareMesh);
  cubes.computeIntersectionAndNormals();
  writeMesh(squareMesh, "results/dualSquares/particles/", 0);
  cubes.extractSurface();

  Carbuncle::TangentParticles2 mergedParticles;
  particles.setDt(0.01745507241); // To rotate 1 degree per timestep
  particles2.setDt(0.01745507241);
  mergedParticles.setDt(0.01745507241);

  mergedParticles =
      Carbuncle::TangentParticles2::mergeParticles(particles, particles2);
  extractLevelsetFromMesh(cubes, squareMesh);
  extractLevelsetFromMesh(cubes2, squareMesh2);
  cubes.merge(cubes2);
  finalMesh = mergedParticles.extractSurface(cubes);

  for (int i = 1; i <= 90; i++) {
    // particles.clearParticles();
    // particles.seedParticlesOverSurface(cubes, squareMesh);
    particles.advectParticles();
    particles2.advectParticles();

    squareMesh = particles.extractSurface(cubes);
    squareMesh2 = particles2.extractSurface(cubes);
    writeMesh(squareMesh, "results/dualSquares/particles/", i);
    writeMesh(squareMesh2, "results/dualSquares/particles2/", i);
    resetLevelset(cubes);
    extractLevelsetFromMesh(cubes, squareMesh);
    extractLevelsetFromMesh(cubes2, squareMesh2);
    cubes.merge(cubes2);
    cubes.computeCellsGradient();
    cubes.computeIntersectionAndNormals();
    cubes.extractSurface();
    printLevelset(cubes);
    cubes.print();

    if (false) {
      Carbuncle::TangentParticles2 mergedParticles;
      mergedParticles =
          Carbuncle::TangentParticles2::mergeParticles(particles, particles2);
      removeParticlesFromSquare(mergedParticles, center, 0.99 * radius, i);
      removeParticlesFromSquare(mergedParticles, -center, 0.99 * radius, i);

      // resampleParticlesFromLevelset(mergedParticles, cubes);

      printParticles(mergedParticles);
      printParticlesPerCell(mergedParticles);
      finalMesh = mergedParticles.extractSurface(cubes);
      writeMesh(finalMesh, "results/dualSquares/analyticRemoval/", i);
    }
    {
      mergedParticles.advectParticles([&] {
        auto &velocity =
            mergedParticles.getParticleArrayData("particleVelocity");
        auto &particleSide =
            mergedParticles.getParticleScalarData("particleSide");
        auto &position =
            mergedParticles.getParticleArrayData("particlePosition");
        for (size_t pid = 0; pid < position.size(); pid++) {
          if (!mergedParticles.isActive(pid))
            continue;
          auto p = position[pid];
          if (particleSide[pid] == 0) {
            p -= center;
          } else {
            p += center;
          }
          velocity[pid] = Eigen::Array2d(-p[1], p[0]);
        }
      });

      resampleParticlesFromLevelset(mergedParticles, cubes);
      resampleParticlesFromGeometry(mergedParticles, cubes, finalMesh);
      printParticles(mergedParticles);
      finalMesh = mergedParticles.extractSurface(cubes);
      writeMesh(finalMesh, "results/dualSquares/fuseParticles/", i);
    }
  }
  // cubes.computeCellsGradient();
  // cubes.computeIntersectionAndNormals();
  // particles.fixLevelsetGradients(cubes);
  // cubes.extractSurface();

  particles.extractSurface(cubes);
  return 0;
}
