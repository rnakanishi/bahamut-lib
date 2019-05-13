#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <graphics/camera.hpp>
#include <graphics/mesh_object.hpp>
#include <gui/event_handler.hpp>
#include <gui/gui.hpp>
#include <iostream>
#include <shader/shader.hpp>

class MyGui : public Garuda::GUI {

public:
  void run() override {

    Garuda::EventHandler events;
    _objects.centerizeObject();
    while (!glfwWindowShouldClose(_window)) {

      // Comandos de entrada
      events.processKeyboardInputs(_window);

      // Comandos de renderizacao vao aqui
      glClearColor(0.7f, 0.75f, 0.75f, 1.0f);
      // glClear(GL_COLOR_BUFFER_BIT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // etc...
      _shader.useShader();

      _objects.draw(_shader);

      // Controla eventos e troca os buffers para renderizacao
      glfwSwapBuffers(_window);
      glfwPollEvents();
    }
  }

  void loadScene() {
    _shader.loadVertexShader("./assets/shaders/camera.vert");
    _shader.loadFragmentShader("./assets/shaders/texture.frag");

    _objects.initialize();
    _objects.loadObjMesh("./assets/3d_models/newdog.obj");
    // _objects.loadObjMesh("./assets/3d_models/dragon.obj");
    // _objects.loadObjMesh("./obj/tetra.obj");
    _objects.loadTexture();
  }

protected:
  Garuda::Shader _shader;
  Garuda::MeshObject _objects;
};

int main(int argc, char const *argv[]) {
  MyGui interface;

  interface.createWindow();
  interface.loadScene();

  interface.run();

  return 0;
}