#include <shader/shader.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace Garuda {
Shader::Shader() {}

Shader::Shader(const char *vertexShaderPath, const char *fragmentShaderPath) {
  vertexShader = loadShader(vertexShaderPath, GL_VERTEX_SHADER);
  fragmentShader = loadShader(fragmentShaderPath, GL_FRAGMENT_SHADER);
  _linkProgram();
}

void Shader::_linkProgram() {
  programId = glCreateProgram();
  glAttachShader(programId, vertexShader);
  glAttachShader(programId, fragmentShader);
  glLinkProgram(programId);

  int success;
  glGetProgramiv(programId, GL_LINK_STATUS, &success);
  if (!success) {
    char infoLog[512];
    glGetProgramInfoLog(programId, 512, NULL, infoLog);
    std::cerr << "Shader::_linkProgram: " << infoLog << std::endl;
  }
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);
}

int Shader::loadShader(const char *shaderPath, unsigned int shaderType) {
  std::string code;
  const char *shaderCode;
  std::ifstream file;
  file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    file.open(shaderPath);
    std::stringstream stream;
    stream << file.rdbuf();
    file.close();
    shaderCode = stream.str().c_str();
  } catch (std::ifstream::failure err) {
    std::cerr << "Shader::loadShader: could not read " << shaderPath
              << std::endl;
  }

  unsigned int shader = -1;
  shader = glCreateShader(shaderType);
  glShaderSource(shader, 1, &shaderCode, NULL);
  glCompileShader(shader);

  int success;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
  if (!success) {
    char infoLog[512];
    glGetShaderInfoLog(shader, 512, NULL, infoLog);
    std::cerr << "Shader::loadShader: shader compilation failed:\n\t" << infoLog
              << std::endl;
  }

  return shader;
}

void Shader::useShader() { glUseProgram(programId); }

int Shader::loadVertexShader(const char *vertexShaderPath) {
  return loadShader(vertexShaderPath, GL_VERTEX_SHADER);
}

int Shader::loadFragmentShader(const char *fragmentShaderPath) {
  return loadShader(fragmentShaderPath, GL_FRAGMENT_SHADER);
}

int Shader::loadGeometryShader(const char *geometryShaderPath) {
  return loadShader(geometryShaderPath, GL_GEOMETRY_SHADER);
}

} // namespace Garuda
