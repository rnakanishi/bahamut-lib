#ifndef __RAMUH_FILE_WRITER_H__
#define __RAMUH_FILE_WRITER_H__

#include <structures/levelset2.h>
#include <string>

namespace Ramuh {

class FileWriter {

public:
  FileWriter();

  void writeLevelSet(LevelSet2 data, const std::string &filename);

  void setDebug(bool out);

private:
  bool _stdout;
};

} // namespace Ramuh

#endif