#version 330 core
out vec4 FragColor;
in vec4 vNormal;
in vec4 fragPosition;

in vec2 texCoord;
uniform int hasTexture;
uniform sampler2D ourTexture;

uniform vec3 ambientLight;
uniform vec3 camPosition;

in vec4 lightPosition;

void main() {
  float ka = 0.5f;
  float kd = 0.5f;
  float ks = 1.0f;
  vec4 lightColor = vec4(1.f);

  // Componente ambiente
  vec4 ambient = vec4(ambientLight, 1.0);

  // componente difuso
  // kd * Il * ( N . L )
  vec4 lightDirection = normalize(lightPosition - fragPosition);
  float diff = max(dot(vNormal, lightDirection), 0.0);
  vec4 diffuse = kd * diff * lightColor;

  // Componente especular
  vec4 viewDirection = normalize(vec4(camPosition, 1.f) - fragPosition);
  vec4 reflectDirection = vec4(reflect(lightDirection, vNormal));
  float spec = pow(max(dot(viewDirection, reflectDirection), 0.0), 64);
  vec4 specular = ks * spec * lightColor;

  if (hasTexture == 1)
    FragColor = (lightColor)*texture(ourTexture, texCoord);
  else
    FragColor = ka * ambient + kd * diffuse + ks * specular;
}