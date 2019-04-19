#ifndef __GARUDA_EVENT_HANDLER_HPP__
#define __GARUDA_EVENT_HANDLER_HPP__
#include <glad/glad.h>
#include <GLFW/glfw3.h>

namespace Garuda {

class EventHandler {
public:
  EventHandler();

  void processKeyboardInputs(GLFWwindow *window);

  void processMouseInputs(GLFWwindow *window);

protected:
  bool _drawPolygons;
};

} // namespace Garuda

#endif