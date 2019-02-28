#ifndef __RAMUH_FILE_WRITER_H__
#define __RAMUH_FILE_WRITER_H__

#include <structures/levelset.h>
#include <string>

namespace Ramuh {

class FileWriter {

public:
  FileWriter();

  static void writeLevelSet(LevelSet data, const std::string &filename);
};

} // namespace Ramuh

#endif