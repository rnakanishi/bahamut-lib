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
  _window = glfwCreateWindow(800, 800, "CG 2019", NULL, NULL);
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
  glViewport(0, 0, 800, 800);

  _shader.loadVertexShader("./assets/vertex_shader.vert");
  _shader.loadFragmentShader("./assets/fragment_shader.frag");

  _objects.initialize();
  _objects.loadObjMesh("./assets/3d_models/dog2.obj");
}

void GUI::run() {

  while (!glfwWindowShouldClose(_window)) {

    // Comandos de entrada
    EventHandler::processKeyboardInputs(_window);

    // Comandos de renderizacao vao aqui
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // etc...
    _objects.draw(_shader);

    // Controla eventos e troca os buffers para renderizacao
    glfwSwapBuffers(_window);
    glfwPollEvents();
  }
}
} // namespace Garuda