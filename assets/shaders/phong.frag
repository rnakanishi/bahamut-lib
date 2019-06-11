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

struct Material {
  vec4 ka, kd, ks;
  float ns;
};

uniform Material material;

void main() {
  vec4 lightColor = vec4(1.f);

  // Componente ambiente
  vec4 ambient = material.ka * vec4(ambientLight, 1.0);

  // componente difuso
  // kd * Il * ( N . L )
  vec4 lightDirection = normalize(lightPosition - fragPosition);
  float diff = max(dot(vNormal, lightDirection), 0.0);
  vec4 diffuse = material.kd * diff * lightColor;

  // Componente especular
  vec4 viewDirection = normalize(vec4(camPosition, 1.f) - fragPosition);
  vec4 reflectDirection = vec4(reflect(-lightDirection, vNormal));
  float spec = pow(max(dot(viewDirection, reflectDirection), 0.0), material.ns);
  vec4 specular = material.ks * spec * lightColor;

  float Kconstant = 1.0;
  float Klinear = 0.07;
  float Kquad = 0.017;
  float lightDistance = length(lightPosition - fragPosition);
  float attenuation = 1.f / (Kconstant + lightDistance * Klinear +
                             (lightDistance * lightDistance) * Kquad);

  if (hasTexture == 1)
    FragColor = (attenuation * (ambient + diffuse + specular)) *
                texture(ourTexture, texCoord);
  else
    FragColor = attenuation * (ambient + diffuse + specular);
  ;
}