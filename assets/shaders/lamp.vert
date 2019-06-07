#version 330 core
layout(location = 0) in vec3 aPos;

// specify a normal output to the fragment shader
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main() {
  gl_Position = projection * view * model * vec4(aPos, 1.f);
  gl_PointSize = 10.0;
}