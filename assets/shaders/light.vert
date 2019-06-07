#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTex;

// specify a normal output to the fragment shader
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 lightTransform;

uniform vec3 light;

out vec2 texCoord;
out vec4 fragPosition;
out vec4 vNormal;
out vec4 lightPosition;

void main() {
  gl_Position = projection * view * model * vec4(aPos, 1.f);
  fragPosition = view * model * vec4(aPos, 1.f);
  vNormal = normalize(vec4(aNormal, 1.f));

  lightPosition = projection * view * lightTransform * vec4(light, 1.f);

  texCoord = aTex;
}