#include <gui/gui.hpp>
#include <gui/event_handler.hpp>

namespace Garuda {
GUI::GUI() {
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
}

GUI::~GUI() { glfwTerminate(); }

void GUI::createWindow() {
  _width = 600;
  _height = 600;
  _window = glfwCreateWindow(_width, _height, "CG 2019", NULL, NULL);
  if (_window == NULL) {
    std::cout << "Failed to create GLFW window\n";
    glfwTerminate();
    exit(-1);
  }
  glfwMakeContextCurrent(_window);
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD\n";
    exit(-2);
  }
  glViewport(0, 0, _width, _height);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_PROGRAM_POINT_SIZE);
}

void GUI::changeViewport(int width, int height) {
  _width = width;
  _height = height;
  glViewport(0, 0, _width, _height);
}

void GUI::run() {

  while (!glfwWindowShouldClose(_window)) {
    // Comandos de renderizacao vao aqui
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Controla eventos e troca os buffers para renderizacao
    glfwSwapBuffers(_window);
    glfwPollEvents();
  }
}
} // namespace Garuda