#version 330 core
out vec4 FragColor;

float near = 0.1;
float far = 10.0;

float backToView(float depth) {
  float z = depth * 2.0 - 1.0;
  // return (z * (far - near) + (far + near)) / 2;
  return (2.0 * near * far) / (far + near - z * (far - near));
}

void main() {
  float depth = backToView(gl_FragCoord.z) / far;
  FragColor = vec4(vec3(depth), 1.0);
}
