#version 330 core
out vec4 FragColor;
in vec4 vNormal;
in vec4 fragPosition;

in vec2 texCoord;
uniform int hasTexture;
uniform sampler2D ourTexture;

uniform vec3 ambientLight;
uniform vec3 camPosition;

#define MAX_LIGHTS 100
uniform mat4 lightTransform[MAX_LIGHTS];
struct LightStruct {
  vec3 position;
  vec4 color;
};
uniform LightStruct lights[MAX_LIGHTS];
uniform int n_activeLights;

struct Material {
  vec4 ka, kd, ks;
  float ns;
};

uniform Material material;

vec4 computeLight(vec4 lightPos, vec4 lightColor) {

  // Componente ambiente
  vec4 ambient = material.ka * vec4(ambientLight, 1.0);

  // componente difuso
  // kd * Il * ( N . L )
  vec4 lightDirection = normalize(lightPos - fragPosition);
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
  float lightDistance = length(lightPos - fragPosition);
  float attenuation = 1.f / (Kconstant + lightDistance * Klinear +
                             (lightDistance * lightDistance) * Kquad);

  return attenuation * (ambient + diffuse + specular);
}

void main() {
  FragColor = vec4(0.f);
  for (int i = 0; i < n_activeLights; i++) {
    vec4 lightPosition = lightTransform[i] * vec4(lights[i].position, 1.f);
    FragColor += computeLight(lightPosition, lights[i].color);
  }

  if (hasTexture == 1)
    FragColor = FragColor * texture(ourTexture, texCoord);
}