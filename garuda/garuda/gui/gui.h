#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>

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

  /// \TODO: create window resizing function

  ///
  /// Run the application in a loop.
  void run();

protected:
  void processInput(GLFWwindow *window);

  int _viewPortSize;
  GLFWwindow *_window;
};
}