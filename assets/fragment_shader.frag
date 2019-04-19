#version 330 core
out vec4 FragColor;
in vec4 vertexColor;
// in vec4 vertexColor; // we set this variable in the OpenGL code.

void main() { FragColor = vec4(0.5, 0.5, 0.5, 1.0); }