#include <utils/mesh_reader.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tiny_obj_loader.h>
#include <glm/vec3.hpp>

namespace Garuda {

void MeshReader::readObj(MeshObject &structure, const char *path) {
  std::string warning, error;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;
  tinyobj::attrib_t attributes;
  bool ret = tinyobj::LoadObj(&attributes, &shapes, &materials, &warning,
                              &error, path);

  for (int i = 0; i < attributes.vertices.size() / 3; i++) {
    structure.addVertex(glm::vec3(attributes.vertices[3 * i + 0],
                                  attributes.vertices[3 * i + 1],
                                  attributes.vertices[3 * i + 2]));
  }
  for (auto shape : shapes)
    for (int i = 0; i < shape.mesh.indices.size() / 3; i++) {
      structure.addFace(glm::ivec3(shape.mesh.indices[3 * i + 0].vertex_index,
                                   shape.mesh.indices[3 * i + 1].vertex_index,
                                   shape.mesh.indices[3 * i + 2].vertex_index));
    }

  std::cout << "Read " << structure.getVerticesSize() << " vertices and "
            << structure.getFacesSize() << " faces\n";
  //   for (auto vertex : attributes.vertices) {
  // structure.addVertex(glm::vec3(vertex[0], vertex[1], vertex[2]));
  //   }
}
} // namespace Garuda