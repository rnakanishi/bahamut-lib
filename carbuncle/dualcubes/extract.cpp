#include <geometry/bounding_box.h>
#include <Eigen/Dense>
#include <utils/file_writer.h>
#include <sstream>
#include <iomanip>
#include "dual_cubes.h"
#include <fluids/levelset_fluid2.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <pugixml.hpp>

class Surfacer3 : public Leviathan::LevelSetFluid3 {
public:
  Surfacer3(Eigen::Array3i gridSize, Ramuh::BoundingBox3 domain)
      : LevelSetFluid3(gridSize, domain) {
    baseFolder = "./";
    nFiles = 1;
    lscount = 0;
    // _dt = 1 / 30.;
  }

  void setBaseFolder(std::string folderName) { baseFolder = folderName; }

  void setNumberOfFiles(int n) { nFiles = n; }

  int getNFiles() { return nFiles; }

  void readFromFile() {
    auto &phi = getCellScalarData(_phiId);
    std::ifstream file;
    std::stringstream filename;
    filename << baseFolder << "ls" << lscount++;
    file.open(filename.str().c_str(), std::ifstream::in);

    if (!file.is_open()) {
      std::cerr << "\033[1;31mError\033[0m Faile to open " << filename.str()
                << std::endl;
      return;
    }

    for (size_t i = 0; i < cellCount(); i++) {
      Eigen::Array3d pos = getCellPosition(i);
      file << pos[0] << " " << pos[1] << " " << pos[2] << " " << phi[i] << "\n";
    }
    std::cout << "File written: " << filename.str() << std::endl;
    file.close();
  }

protected:
  std::string baseFolder;
  int lscount, nFiles;
};

class Configuration {
public:
  Configuration();

  void configureFromFile(std::string filename, Surfacer3 &surface) {
    auto parseResult = _document.load_file(filename.c_str());
    auto root = _document.children("simulation");

    _numberOfSimulations = std::distance(root.begin(), root.end());
    _currentSimulation = 0;

    pugi::xml_node node;
    node = _document.child("simulation");
    for (int siblings = 0; siblings < _currentSimulation; siblings++)
      node = node.next_sibling();

    node = node.child("parameters");

    // Simulatino parameters
    int size, frames;
    double domainMin, domainMax, dt = 1 / 60.;
    int maxParticles;

    sscanf(node.child("gridSize").child_value(), "%d", &size);
    sscanf(node.child("domain").child_value(), "%lf %lf", &domainMin,
           &domainMax);
    sscanf(node.child("frames").child_value(), "%d", &frames);

    Eigen::Array3i gridSize(size);
    Ramuh::BoundingBox3 domain(domainMin, domainMax);

    surface.setNumberOfFiles(frames);
    surface.setGridSize(Eigen::Array3i(size));
    surface.setDomain(Ramuh::BoundingBox3(domainMin, domainMax));
  }

protected:
  int _numberOfSimulations;
  int _currentSimulation;
  pugi::xml_document _document;
};

int main(int argc, char const *argv[]) {

  Surfacer3 surface(Eigen::Array3i(128), Ramuh::BoundingBox3(0, 1));

  Configuration config;
  config.configureFromFile("./carbuncle/configs/pls.xml", surface);

  for (size_t i = 0; i < surface.getNFiles(); i++) {
    /* code */
    surface.readFromFile();
    surface.extractSurface();
  }

  Carbuncle::DualCubes3 cubes(
      Eigen::Array3i(40, 40, 40),
      Ramuh::BoundingBox3(Eigen::Array3d(-5, -5, -5), Eigen::Array3d(5, 5, 5)));
  cubes.initialize(Eigen::Array3d(0, 0, 0), 3,
                   Carbuncle::DualCubes3::ParametricSurface::ELLIPSOID);
  cubes.defineVelocity();
  // cubes.printCells();

  // cubes.computeIntersection();
  // cubes.analyticNormals(Eigen::Array3d(0, 0, 0), 0.5,
  //                       Carbuncle::DualCubes3::ParametricSurface::CUBE);
  // cubes.extractSurface();
  // return 1;

  // cubes.computeCellsGradient();
  // cubes.computeIntersectionAndNormals();
  // cubes.extractSurface();
  cubes.redistance();
  cubes.print();

  return 1;

  for (int i = 1; i <= 150; i++) {
    cubes.advectWeno();
    // cubes.advectUpwind();

    if (i % 10 == 0)
      cubes.redistance();
    // cubes.computeIntersection();
    // cubes.computeNormals();

    cubes.computeCellsGradient();
    cubes.computeIntersectionAndNormals();
    cubes.extractSurface();
    cubes.print();
  }

  //   auto surface = cubes.marc hingTetrahedra();
  //   Ramuh::FileWriter writer;
  //   std::ostringstream objname;
  //   objname << "results/marching/tetra_" << std::setfill('0') << std::setw(4)
  //   << 0
  //   << ".obj";
  //   writer.writeMeshModel(surface, objname.str());
  // }
  return 0;
}
