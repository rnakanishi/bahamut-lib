#include <utils/callbacks.hpp>
#include <graphics/camera.hpp>
#include <iostream>
#include <glm/gtx/string_cast.hpp>

namespace Garuda {

void __mouseButton(GLFWwindow *window, int button, int action, int mods) {
  Camera *camera;
  camera = static_cast<Camera *>(glfwGetWindowUserPointer(window));

  static bool leftPressed = false;
  static bool firstMouse = true;
  static double lastX, lastY;
  switch (button) {
  case GLFW_MOUSE_BUTTON_1:
    if (action == GLFW_PRESS) {
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
      leftPressed = true;
      firstMouse = true;
    } else if (action == GLFW_RELEASE) {
      std::cerr << "Release\n";
      leftPressed = false;
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
    break;
  case GLFW_MOUSE_BUTTON_2:
    std::cerr << "Button 2\n";
    break;
  default:
    std::cerr << "Other button";
  }

  if (leftPressed) {
    double x, y;
    glfwGetCursorPos(window, &x, &y);

    if (firstMouse) {
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
      firstMouse = false;
      lastX = x;
      lastY = y;
    }

    float sensitivity = 0.001f;
    float dx = (x - lastX) * sensitivity;
    float dy = (y - lastY) * sensitivity;
    lastX = x;
    lastY = y;
    std::cerr << x << ' ' << y << std::endl;
    std::cerr << dx << ' ' << dy << std::endl;
    camera->rotateCamera(-dy, -dx);
  }
}

void __mouseCursor(GLFWwindow *window, double x, double y) {
  Camera *camera;
  camera = static_cast<Camera *>(glfwGetWindowUserPointer(window));
  static bool firstMouse = true;
  static double lastX, lastY;
  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
    // double x, y;
    // glfwGetCursorPos(window, &x, &y);
    if (firstMouse) {
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
      firstMouse = false;
      lastX = x;
      lastY = y;
    }

    float sensitivity = 0.001f;
    float dx = (x - lastX) * sensitivity;
    float dy = (y - lastY) * sensitivity;
    lastX = x;
    lastY = y;
    std::cerr << x << ' ' << y << std::endl;
    std::cerr << dx << ' ' << dy << std::endl;
    camera->rotateCamera(-dy, -dx);
  }
  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE) {
    if (!firstMouse)
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    firstMouse = true;
  }
}

void __mouseScroll(GLFWwindow *window, double dx, double dy) {}

} // namespace Garuda