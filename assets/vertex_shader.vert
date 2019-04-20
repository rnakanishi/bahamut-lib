#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 color;

out vec4 vertexColor; // specify a color output to the fragment shader
uniform mat4 transform;

void main() { gl_Position = transform * vec4(aPos, 1.0f); }