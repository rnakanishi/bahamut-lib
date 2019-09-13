#ifndef __RAMUH_FILE_WRITER_H__
#define __RAMUH_FILE_WRITER_H__

#include <structures/triangle_mesh.h>
#include <string>

namespace Ramuh {

class FileWriter {

public:
  FileWriter();

  void writeMeshModel(TriangleMesh model, const std::string &filename);

  void setDebug(bool out);

private:
  bool _stdout;
};

} // namespace Ramuh

#endif