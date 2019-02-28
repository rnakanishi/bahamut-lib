#include <utils/file_writer.h>
#include <fstream>

namespace Ramuh {
FileWriter::FileWriter() {}

void FileWriter::writeLevelSet(LevelSet data, const std::string &filename) {
  std::ofstream file(filename, std::ofstream::out);

  for (int k = 0; k < data.resolution().x(); k++) {
    for (int j = 0; j < data.resolution().y(); j++) {
      for (int i = 0; i < data.resolution().z(); i++) {
        file << data[i][j][k] << " ";
      }
    }
  }

  file.close();
}

} // namespace Ramuh
