#version 330 core
out vec4 FragColor;

in vec2 texCoord;
uniform int hasTexture;
uniform sampler2D ourTexture;

in vec4 vColor;

void main() {
  vec4 lightColor = vec4(1.f);

  if (hasTexture == 1)
    FragColor = (lightColor)*texture(ourTexture, texCoord);
  else
    FragColor = vColor;
}