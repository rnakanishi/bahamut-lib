#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <shader/shader.hpp>
#include <graphics/mesh_object.hpp>

namespace Garuda {

class GUI {
public:
  ///
  /// Instantiates glfw environment
  GUI();

  ///
  /// Finalizes glfw environment
  ~GUI();

  ///
  /// Create single window of 800x600 size
  /// \TODO: Fix window size to be class attribute
  void createWindow();

  void changeViewport(int width, int height);

  /// \TODO: create window resizing function

  ///
  /// Run the application in a loop.
  virtual void run();

protected:
  int _viewPortSize;
  int _width, _height;
  GLFWwindow *_window;
};
} // namespace Garuda