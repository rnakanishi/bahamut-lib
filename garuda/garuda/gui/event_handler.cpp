#include <gui/event_handler.hpp>

namespace Garuda {
void EventHandler::processKeyboardInputs(GLFWwindow *window) {
  // void EventHandler::processKeyboardInputs(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, true);
}

void EventHandler::processMouseInputs(GLFWwindow *window) {}

} // namespace Garuda