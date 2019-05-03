#include <structures/levelset3.h>
#include <fstream>
#include <iostream>
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
          file >> _phi[i][j][k];
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

  surface.readData("data/26");
  auto triangles = surface.marchingTetrahedra();

  Ramuh::FileWriter writer;
  writer.writeMeshModel(triangles, "obj/reconstruct.obj");
}