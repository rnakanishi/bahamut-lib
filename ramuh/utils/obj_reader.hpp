#ifndef __RAMUH_OBJ_READER_HPP__
#define __RAMUH_OBJ_READER_HPP__

#include <structures/triangle_mesh.h>
#include <tiny_obj_loader.h>

namespace Ramuh {
class ObjReader {
public:
  Ramuh::TriangleMesh readObj(const char *objPath);

protected:
  ObjReader();
};

} // namespace Ramuh

#endif