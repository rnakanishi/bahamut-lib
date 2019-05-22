#ifndef __GARUDA_MESH_READER_HPP__
#define __GARUDA_MESH_READER_HPP__

#include <graphics/mesh_object.hpp>

namespace Garuda {

class MeshReader {
public:
  void readVertices(void *data, float x, float y, float z, float w);

  static void readObj(MeshObject &structure, const char *path);

  static void readObjWithTexture(MeshObject &structure, const char *path);

  static bool readTexture(MeshObject &structure, const char *path);

  static void objLoader(MeshObject &structure, const char *path);
};

} // namespace Garuda

#endif