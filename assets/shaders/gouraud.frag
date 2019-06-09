#version 330 core
out vec4 FragColor;
in vec4 vNormal;
in vec4 fragPosition;

in vec2 texCoord;
uniform int hasTexture;
uniform sampler2D ourTexture;
// in vec4 vertexColor; // we set this variable in the OpenGL code.
in vec4 diffuse;
in vec4 specular;
in vec4 vColor;

void main() {
  float ka = 0.5f;
  float kd = 0.2f;
  float ks = 1.0f;
  vec4 lightColor = vec4(1.f);

  if (hasTexture == 1)
    FragColor = (lightColor)*texture(ourTexture, texCoord);
  else
    FragColor = vColor;
}