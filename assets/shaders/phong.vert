#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTex;

// specify a normal output to the fragment shader
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform vec3 light;

out vec2 texCoord;
out vec4 diffuse;
out vec4 specular;

out vec4 vNormal;
out vec4 fragPosition;

void main() {
  vec4 lightColor = vec4(1.f);

  gl_Position = projection * view * model * vec4(aPos, 1.f);

  fragPosition = model * vec4(aPos, 1.f);
  vec3 normal = mat3(transpose(inverse(model))) * aNormal;
  vNormal = normalize(vec4(normal, 1.f));

  texCoord = aTex;
}