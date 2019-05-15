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

class MyGui : public Garuda::GUI {

public:
  void run() override {
    Garuda::EventHandler events;
    while (!glfwWindowShouldClose(_window)) {

      // Comandos de entrada
      Garuda::EventHandler::processKeyboardInputs(_window);

      // Comandos de renderizacao vao aqui
      glClearColor(0.7f, 0.75f, 0.75f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // etc...
      _scene.draw();

      // Controla eventos e troca os buffers para renderizacao
      glfwSwapBuffers(_window);
      glfwPollEvents();
    }
  }

  void loadScene() { _scene.load(); }

protected:
  Garuda::Scene _scene;
};

int main(int argc, char const *argv[]) {
  MyGui interface;

  interface.createWindow();
  interface.loadScene();

  interface.run();

  return 0;
}