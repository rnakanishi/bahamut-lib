#include <utils/file_writer.h>
#include <glm/vec3.hpp>
#include <fstream>
#include <iostream>

namespace Ramuh {
FileWriter::FileWriter() { _stdout = false; }

void FileWriter::setDebug(bool out) { _stdout = out; }

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

void FileWriter::writeArrayToFile(std::string filename,
                                  std::vector<int> &array) {
  std::ofstream file;
  file.open(filename, std::ofstream::out);
  if (!file.is_open()) {
    std::cerr << "\033[1;31mError\033[0m Timer::logToFile: Failed to open "
              << filename << std::endl;
    return;
  }
  for (auto value : array) {
    file << value << " ";
  }
  file.close();
  std::cout << "File written: " << filename << std::endl;
}

} // namespace Ramuh
