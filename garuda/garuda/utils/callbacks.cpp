#include <graphics/camera.hpp>
#include <utils/callbacks.hpp>
#include <iostream>
#include <glm/gtx/string_cast.hpp>

namespace Garuda {

void __mouseButton(GLFWwindow *window, int button, int action, int mods) {}

void __mouseCursor(GLFWwindow *window, double xPos, double yPos) {
  Camera *camera;
  camera = static_cast<Camera *>(glfwGetWindowUserPointer(window));
  static bool firstMouse = true;
  static double lastX = xPos, lastY = yPos;

  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
    // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    if (firstMouse) {
      firstMouse = false;
      lastX = x;
      lastY = y;
    }

    float sensitivity = 0.001f;
    float dx = (x - lastX) * sensitivity;
    float dy = (y - lastY) * sensitivity;
    lastX = x;
    lastY = y;
    std::cerr << x << " " << y << " " << dx << " " << dy << std::endl;
    camera->rotateCamera(dy, dx);
  }
  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    std::cerr << "Right button pressed\n";
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    if (firstMouse) {
      firstMouse = false;
      lastX = x;
      lastY = y;
    }

    float sensitivity = 1.f;
    float dx = (x - lastX) * sensitivity;
    float dy = (y - lastY) * sensitivity;
    std::cerr << x << " " << y << " " << dx << " " << dy << std::endl;
    lastX = x;
    lastY = y;
    camera->moveCamera(glm::vec3(dx, dy, 0.f));
  } else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) ==
             GLFW_PRESS) {
    std::cerr << "Middle button pressed\n";
  }

  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE &&
      glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE) {
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    firstMouse = true;
  }
}

void __mouseScroll(GLFWwindow *window, double dx, double dy) {
  Camera *camera;
  camera = static_cast<Camera *>(glfwGetWindowUserPointer(window));
  camera->setFrustum((float)dy);
  camera->perspectiveProjection();
}

void __viewportChange(GLFWwindow *window, int width, int height) {
  Camera *camera;
  camera = static_cast<Camera *>(glfwGetWindowUserPointer(window));
  camera->setAspect(width, height);
}

} // namespace Garuda