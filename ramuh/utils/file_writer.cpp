#include <utils/file_writer.h>
#include <fstream>

namespace Ramuh {
FileWriter::FileWriter() { _stdout = false; }

void FileWriter::setDebug(bool out) { _stdout = out; }

void FileWriter::writeLevelSet(LevelSet2 data, const std::string &filename) {
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

} // namespace Ramuh
