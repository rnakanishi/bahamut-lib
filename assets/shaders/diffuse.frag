#version 330 core
out vec4 FragColor;
in vec4 vertexNormal;
in vec4 fragPos;
// in vec4 vertexColor; // we set this variable in the OpenGL code.

void main() {
  vec4 lightPosition = vec4(0.0, 2.0, 2.0, 1.0);
  vec4 normal = normalize(vertexNormal);
  vec4 lightDirection = normalize(lightPosition - fragPos);
  float diff = max(dot(normal, lightDirection), 0.0);

  FragColor = diff * vec4(0.75, 0.75, 0.75, 1.0);
}