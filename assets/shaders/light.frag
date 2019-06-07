#version 330 core
out vec4 FragColor;
in vec4 vNormal;
in vec4 fragPosition;

in vec2 texCoord;
in vec4 lightPosition;
uniform int hasTexture;
uniform sampler2D ourTexture;
// in vec4 vertexColor; // we set this variable in the OpenGL code.
uniform vec3 ambientLight;

void main() {
  vec4 lightDirection = normalize(lightPosition - fragPosition);
  float diff = max(dot(vNormal, lightDirection), 0.0);

  vec4 ambient = vec4(ambientLight, 1.0);

  if (hasTexture == 1)
    FragColor = diff * texture(ourTexture, texCoord);
  else
    FragColor = diff * vec4(1.f);
}