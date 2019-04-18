#version 330 core
layout(location = 0) in vec3 aPos;

out vec4 vertexColor; // specify a color output to the fragment shader

void main() {
  gl_Position = vec4(aPos, 1.0);
  vertexColor = vec4(0.5, 0.151, 0.55, 1.0);
}