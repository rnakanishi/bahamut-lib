#ifndef __GARUDA_EVENT_HANDLER_HPP__
#define __GARUDA_EVENT_HANDLER_HPP__
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <graphics/camera.hpp>

namespace Garuda {

class EventHandler {
public:
  static void processKeyboardInputs(GLFWwindow *window);

  static void cameraKeyboardInputs(GLFWwindow *window, Camera &camera);

  static void cameraMouseInputs(GLFWwindow *window, Camera &camera);

  static void processMouseInputs(GLFWwindow *window);
};

} // namespace Garuda

#endif