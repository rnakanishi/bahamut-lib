#include <utils/file_writer.h>
#include <glm/vec3.hpp>
#include <fstream>

namespace Ramuh {
FileWriter::FileWriter() { _stdout = false; }

void FileWriter::setDebug(bool out) { _stdout = out; }

void FileWriter::writeLevelSet(LevelSet2 &data, const std::string &filename) {
  std::ofstream file;
  file.open(filename, std::ofstream::out);

  if (!file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return;
  }

  std::cout << "Writing to " << filename << std::endl;
  for (int j = 0; j < data.resolution().y(); j++) {
    for (int i = 0; i < data.resolution().x(); i++) {
      file << data[i][j] << " ";
      if (_stdout)
        std::cout << data[i][j] << " ";
    }
    file << std::endl;
    if (_stdout)
      std::cout << std::endl;
  }

  file.close();
}

void FileWriter::writeLevelSet(LevelSet3 &data, const std::string &filename) {
  std::ofstream file;
  file.open(filename, std::ofstream::out);

  if (!file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return;
  }

  std::cout << "Writing to " << filename << std::endl;
  for (int k = 0; k < data.resolution().z(); k++) {
    for (int i = 0; i < data.resolution().x(); i++) {
      for (int j = 0; j < data.resolution().y(); j++) {
        file << data[i][j][k] << " ";
        if (_stdout)
          std::cout << data[i][j][k] << " ";
      }
      // file << std::endl;
      if (_stdout)
        std::cout << std::endl;
    }
  }
  file.close();
}
void FileWriter::writeMeshModel(TriangleMesh model,
                                const std::string &filename) {
  std::ofstream file;
  int nvertices, nfaces;
  nvertices = model.getVerticesSize();
  nfaces = model.getFacesSize();

  file.open(filename, std::ofstream::out);
  if (!file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return;
  }

  std::cout << "Writing model " << filename << std::endl;
  std::cout << nvertices << " vertices and " << nfaces << " faces\n";

  glm::vec3 vertex;
  glm::ivec3 face;
  for (int i = 0; i < nvertices; i++) {
    vertex = model.getVertex(i);
    file << "v ";
    file << vertex[0] << ' ' << vertex[1] << ' ' << vertex[2];
    file << std::endl;
  }
  file << std::endl;
  for (int i = 0; i < nfaces; i++) {
    face = model.getFace(i);
    file << "f ";
    file << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1;
    file << std::endl;
  }
  file.close();
}

} // namespace Ramuh
