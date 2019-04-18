#ifndef __GARUDA_EVENT_HANDLER_HPP__
#define __GARUDA_EVENT_HANDLER_HPP__
#include <glad/glad.h>
#include <GLFW/glfw3.h>

namespace Garuda {

class EventHandler {
public:
  static void processKeyboardInputs(GLFWwindow *window);

  static void processMouseInputs(GLFWwindow *window);
};
void EventHandler::processKeyboardInputs(GLFWwindow *window) {
  // void EventHandler::processKeyboardInputs(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, true);
}

} // namespace Garuda

#endif