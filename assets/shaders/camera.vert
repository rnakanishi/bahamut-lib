#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec2 aTex;

// specify a normal output to the fragment shader
uniform mat4 model;
out vec2 texCoord;

void main() {
  gl_Position = model * vec4(aPos, 1.0f);
  texCoord = aTex;
}