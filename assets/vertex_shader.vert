#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 normal;

// specify a normal output to the fragment shader
out vec4 vertexNormal, fragPos;
uniform mat4 transform;

void main() {
  gl_Position = transform * vec4(aPos, 1.0f);
  fragPos = gl_Position;
  vertexNormal = vec4(normal, 1.0);
}