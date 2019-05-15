#include <gui/event_handler.hpp>

namespace Garuda {
EventHandler::EventHandler() { _drawPolygons = false; }

void EventHandler::processKeyboardInputs(GLFWwindow *window) {
  // void EventHandler::processKeyboardInputs(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, true);
  if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void EventHandler::processMouseInputs(GLFWwindow *window) {}

void EventHandler::cameraKeyboardInputs(Camera &camera) {}

void EventHandler::cameraMouseInputs(Camera &camera) {}

} // namespace Garuda