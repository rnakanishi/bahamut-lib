#include <iostream>
#include <gui/gui.hpp>
#include <gui/event_handler.hpp>
#include <shader/shader.hpp>
#include <graphics/mesh_object.hpp>
#include <structures/levelset3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class LevelSetObject : public Garuda::MeshObject {
public:
  void loadMesh() {

    glBindVertexArray(_vao);
    // Assign vertex position buffer
    glBindBuffer(GL_ARRAY_BUFFER, _vbo[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * getVerticesSize(),
                 _vertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3),
                          (void *)0);
    glEnableVertexAttribArray(0);

    // Assign vertices indices to triangles
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(glm::ivec3) * getFacesSize(),
                 _faces.data(), GL_STATIC_DRAW);
    std::cout << "Vertice:";
    std::cout << "Vertice:" << _vertices.size();
    // // Assign normal to each vertex
    // glBindBuffer(GL_ARRAY_BUFFER, _vbo[1]);
    // glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * getVerticesSize(),
    //              _vertexNormal.data(), GL_STATIC_DRAW);
    // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3),
    //                       (void *)0);
    glEnableVertexAttribArray(1);

    // Ensure no buffer is binded
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
  }

  void draw(Garuda::Shader shader) override {
    shader.useShader();
    glBindVertexArray(_vao);
    // glDrawArrays(GL_TRIANGLES, 0, 3 * getVerticesSize());
    glDrawElements(GL_TRIANGLES, 3 * getFacesSize(), GL_UNSIGNED_INT, 0);
  }

protected:
  Ramuh::LevelSet3 _levelset;
};

class MyGui : public Garuda::GUI {

public:
  void run() override {

    Garuda::EventHandler events;
    glm::vec3 size = _objects->getBBoxSize(),
              centroid = _objects->getBboxCenter();
    // double largerDimension = std::max(std::max(size[0], size[1]), size[2]);
    // std::cout << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
    // std::cout << centroid[0] << ' ' << centroid[1] << ' ' << centroid[2]
    //           << std::endl;
    // std::cout << largerDimension << std::endl;

    while (!glfwWindowShouldClose(_window)) {

      // Comandos de entrada
      events.processKeyboardInputs(_window);

      // Comandos de renderizacao vao aqui
      glClearColor(0.7f, 0.75f, 0.75f, 1.0f);
      // glClear(GL_COLOR_BUFFER_BIT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // etc...
      _shader.useShader();
      // glm::mat4 trans = glm::mat4(1.0f);
      // glm::translate(trans, centroid);
      // trans =
      //     glm::rotate(trans, (float)glfwGetTime(), glm::vec3(0.0f, 1.0f,
      //     0.0f));
      // trans = glm::scale(trans, glm::vec3(1.5 / largerDimension));
      // trans = glm::translate(trans, -centroid);

      // unsigned int transformLoc =
      //     glGetUniformLocation(_shader.getId(), "transform");
      // glUniformMatrix4fv(transformLoc, 1, GL_FALSE, glm::value_ptr(trans));

      _objects->draw(_shader);

      // Controla eventos e troca os buffers para renderizacao
      glfwSwapBuffers(_window);
      glfwPollEvents();
    }
  }

  void loadScene() {
    _shader.loadVertexShader("./assets/shaders/texture.vert");
    _shader.loadFragmentShader("./assets/shaders/texture.frag");

    Ramuh::LevelSet3 _levelset;
    int resolution = 24;
    _levelset.setResolution(Ramuh::Vector3i(resolution));
    _levelset.setSize(Ramuh::Vector3d(1.0, 1.0, 1.0));
    _levelset.addSphereSurface(Ramuh::Vector3d(0.5, 0.5, 0.5), 0.25);

    Ramuh::TriangleMesh surface = _levelset.marchingTetrahedra();
    _objects = (LevelSetObject *)&surface;

    _objects->loadMesh();
  }

protected:
  Garuda::Shader _shader;
  LevelSetObject *_objects;
};

int main(int argc, char const *argv[]) {
  MyGui interface;

  interface.createWindow();

  interface.loadScene();
  interface.run();

  return 0;
}