#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTex;

// specify a normal output to the fragment shader
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform vec3 camPosition;

uniform vec3 light;
uniform mat4 lightTransform;
uniform vec3 ambientLight;

out vec2 texCoord;
out vec4 diffuse;
out vec4 specular;

out vec4 vColor;

void main() {
  texCoord = aTex;

  float ka = 0.5;
  float kd = 0.3;
  float ks = 1.0;

  gl_Position = projection * view * model * vec4(aPos, 1.f);

  // Componente ambiente
  vec4 ambient = ka * vec4(ambientLight, 1.0);
  vec4 lightColor = vec4(1.f);

  vec4 lightPosition = lightTransform * vec4(light, 1.f);
  vec4 vPosition = model * vec4(aPos, 1.f);

  // Componente difuso
  vec4 vNormal = normalize(vec4(aNormal, 1.f));
  vec4 lightDirection = normalize(lightPosition - vPosition);
  float diff = max(dot(vNormal, lightDirection), 0.0);
  diffuse = kd * diff * lightColor;

  // Componente especular
  vec4 viewDirection = normalize(vec4(camPosition, 1.f) - vPosition);
  vec4 reflectDirection = reflect(-lightDirection, vNormal);
  float spec = pow(max(dot(viewDirection, reflectDirection), 0.0), 32);
  specular = ks * spec * lightColor;

  float Kconstant = 1.0;
  float Klinear = 0.09;
  float Kquad = 0.032;
  float lightDistance = length(lightPosition - vPosition);
  float attenuation = 1.f / (Kconstant + lightDistance * Klinear +
                             (lightDistance * lightDistance) * Kquad);

  vColor = attenuation * (ambient + diffuse + specular);
}