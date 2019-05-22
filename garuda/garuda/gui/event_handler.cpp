#include <gui/event_handler.hpp>
#include <iostream>

namespace Garuda {

void EventHandler::processKeyboardInputs(GLFWwindow *window) {
  // void EventHandler::processKeyboardInputs(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, true);
  if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void EventHandler::processMouseInputs(GLFWwindow *window) {
  glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
  bool firstMouse = true;
  int lastX, lastY;
  float yaw = -90.f;
  float pitch = 0.f;
  auto mouseCallback = [&](GLFWwindow *window, double xpos, double ypos) {
    if (firstMouse) {
      firstMouse = false;
      lastX = xpos;
      lastY = ypos;
    }
    float sensitivity = 0.01f;
    float dx = (xpos - lastX) * sensitivity;
    float dy = (ypos - lastY) * sensitivity;
  };

  // glfwSetCursorPosCallback(window, mouseCallback);
}

void EventHandler::cameraKeyboardInputs(GLFWwindow *window, Camera &camera) {
  if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
    camera.moveCamera(glm::vec3(0.f, 0.f, -1.f));
  if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
    camera.moveCamera(glm::vec3(0.f, 0.f, 1.f));
  if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
    camera.moveCamera(glm::vec3(-1.f, 0.f, 0.f));
  if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
    camera.moveCamera(glm::vec3(1.f, 0.f, 0.f));
  if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
    camera.moveCamera(glm::vec3(0.f, 1.f, 0.f));
  if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
    camera.moveCamera(glm::vec3(0.f, -1.f, 0.f));

  if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
    camera.orthogonalProjection();
  if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
    camera.perspectiveProjection();

  static bool firstMouse = true;
  static double lastX, lastY;
  if (glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS) {
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

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
    std::cerr << x << ' ' << y << std::endl;
    std::cerr << dx << ' ' << dy << std::endl;
    camera.rotateCamera(dy, dx);
  }
  if (glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_RELEASE) {
    firstMouse = true;
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
  }
}

void EventHandler::cameraMouseInputs(GLFWwindow *window, Camera &camera) {
  static bool firstMouse = true;
  static double lastX, lastY;
  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
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
    camera.rotateCamera(dy, dx);
  }
  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE) {
    if (!firstMouse)
      glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    firstMouse = true;
  }
}

} // namespace Garuda