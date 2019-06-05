#include <structures/levelset3.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utils/file_writer.h>

class LevelSetTriangles : public Ramuh::LevelSet3 {
public:
  LevelSetTriangles() {}

  void readData(const char *filename) {
    std::ifstream file;
    file.open(filename);
    if (!file.is_open()) {
      std::cerr << "Error opening " << filename << std::endl;
      return;
    }
    for (int k = 0; k < _resolution[2]; k++) {
      for (int i = 0; i < _resolution[0]; i++) {
        for (int j = 0; j < _resolution[1]; j++) {
          file >> _phi[_currBuffer][i][j][k];
        }
      }
    }
  }

protected:
};

int main(void) {
  LevelSetTriangles surface;

  int resolution = 64;
  surface.setResolution(Ramuh::Vector3i(resolution));
  surface.setSize(Ramuh::Vector3d(1.0, 1.0, 1.0));

  for (int frame = 0; frame < 1; frame++) {
    std::stringstream objname;
    std::stringstream dataname;
    objname << "results/marching/" << std::setw(4) << std::setfill('0') << frame
            << ".obj";
    dataname << "results/datatest/" << frame;

    surface.readData(dataname.str().c_str());
    auto triangles = surface.marchingTetrahedra();

    Ramuh::FileWriter writer;
    writer.writeMeshModel(triangles, objname.str().c_str());
  }
}