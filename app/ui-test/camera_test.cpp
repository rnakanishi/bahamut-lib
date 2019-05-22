#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <graphics/camera.hpp>
#include <graphics/mesh_object.hpp>
#include <graphics/scene.hpp>
#include <gui/event_handler.hpp>
#include <gui/gui.hpp>
#include <iostream>
#include <shader/shader.hpp>
#include <GLFW/glfw3.h>
#include <utils/callbacks.hpp>
#include <glm/gtx/string_cast.hpp>

class MyGui : public Garuda::GUI {

public:
  void run() override {
    Garuda::EventHandler events;
    _scene.getActiveCamera().setLookAt(glm::vec3(0.0f));
    while (!glfwWindowShouldClose(_window)) {

      // Comandos de entrada
      Garuda::EventHandler::processKeyboardInputs(_window);
      Garuda::EventHandler::cameraKeyboardInputs(_window,
                                                 _scene.getActiveCamera());
      // Garuda::EventHandler::cameraMouseInputs(_window,
      // _scene.getActiveCamera());

      // Comandos de renderizacao vao aqui
      glClearColor(1.f, 1.f, 1.f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // etc...
      _scene.draw();

      // Controla eventos e troca os buffers para renderizacao
      glfwSwapBuffers(_window);
      glfwPollEvents();
    }
  }

  void loadScene() {
    _scene.getActiveCamera().setAspect(600, 600);
    _scene.load();
  }

  void prepareCallbacks() {
    glfwSetWindowUserPointer(_window, &_scene.getActiveCamera());
    // glfwSetMouseButtonCallback(_window, Garuda::__mouseButton);
    glfwSetCursorPosCallback(_window, Garuda::__mouseCursor);
    glfwSetScrollCallback(_window, Garuda::__mouseScroll);
    glfwSetFramebufferSizeCallback(_window, Garuda::__viewportChange);
  }

protected:
  Garuda::Scene _scene;
};

int main(int argc, char const *argv[]) {
  MyGui interface;

  interface.createWindow();
  interface.prepareCallbacks();
  interface.loadScene();

  interface.run();

  return 0;
}