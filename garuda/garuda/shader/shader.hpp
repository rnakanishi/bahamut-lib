#ifndef __GARUDA_SHADER_HPP__
#define __GARUDA_SHADER_HPP__
#include <glad/glad.h>

namespace Garuda {
class Shader {

public:
  Shader();
  Shader(const char *vertexShaderPath, const char *fragmentShaderPath);

  /**
   * @brief Use current program shader
   *
   **/
  void useShader();

  int loadShader(const char *shaderPath, unsigned int shaderType);

  int loadVertexShader(const char *geometryShaderPath);

  int loadGeometryShader(const char *geometryShaderPath);

  int loadFragmentShader(const char *fragmentShaderPath);

protected:
  void _linkProgram();

  unsigned int programId, vertexShader, fragmentShader;
  bool needLink;
};

} // namespace Garuda

#endif