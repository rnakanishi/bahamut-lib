#include <iostream>
#include <gui/gui.hpp>
#include <gui/event_handler.hpp>
#include <shader/shader.hpp>
#include <graphics/mesh_object.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class MyGui : public Garuda::GUI {

public:
  void run() override {

    Garuda::EventHandler events;
    glm::vec3 size = _objects.getBBoxSize(),
              centroid = _objects.getBboxCenter();
    double largerDimension = std::max(std::max(size[0], size[1]), size[2]);
    std::cout << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
    std::cout << centroid[0] << ' ' << centroid[1] << ' ' << centroid[2]
              << std::endl;
    std::cout << largerDimension << std::endl;

    while (!glfwWindowShouldClose(_window)) {

      // Comandos de entrada
      events.processKeyboardInputs(_window);

      // Comandos de renderizacao vao aqui
      glClearColor(0.7f, 0.75f, 0.75f, 1.0f);
      // glClear(GL_COLOR_BUFFER_BIT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // etc...
      _shader.useShader();
      glm::mat4 trans = glm::mat4(1.0f);
      glm::translate(trans, centroid);
      trans =
          glm::rotate(trans, (float)glfwGetTime(), glm::vec3(0.0f, 1.0f, 0.0f));
      trans = glm::scale(trans, glm::vec3(1.5 / largerDimension));
      trans = glm::translate(trans, -centroid);

      unsigned int transformLoc =
          glGetUniformLocation(_shader.getId(), "transform");
      glUniformMatrix4fv(transformLoc, 1, GL_FALSE, glm::value_ptr(trans));

      _objects.draw(_shader);

      // Controla eventos e troca os buffers para renderizacao
      glfwSwapBuffers(_window);
      glfwPollEvents();
    }
  }

  void loadScene() {
    _shader.loadVertexShader("./assets/shaders/texture.vert");
    _shader.loadFragmentShader("./assets/shaders/texture.frag");

    _objects.initialize();
    // _objects.loadTexture();
    _objects.loadObjMesh("./assets/3d_models/dragonSimple.obj");
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