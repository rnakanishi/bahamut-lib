#version 330 core
out vec4 FragColor;

in vec2 texCoord;
uniform sampler2D ourTexture;
// in vec4 vertexColor; // we set this variable in the OpenGL code.

void main() { FragColor = texture(ourTexture, texCoord); }