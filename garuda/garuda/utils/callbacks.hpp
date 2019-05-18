#ifndef __GARUDA_UTILS_CALLBACKS_HPP
#define __GARUDA_UTILS_CALLBACKS_HPP
#include <GLFW/glfw3.h>

namespace Garuda {
/**
 * @brief
 *
 * @param window
 * @param button
 * @param action
 * @param mods
 */
void __mouseButton(GLFWwindow *window, int button, int action, int mods);

/**
 * @brief
 *
 * @param window
 * @param x
 * @param y
 */
void __mouseCursor(GLFWwindow *window, double x, double y);

/**
 * @brief
 *
 * @param window
 * @param dx
 * @param dy
 */
void __mouseScroll(GLFWwindow *window, double dx, double dy);

/**
 * @brief
 *
 * @param window
 * @param width
 * @param height
 **/
void __viewportChange(GLFWwindow *window, int width, int height);

} // namespace Garuda
#endif // !__GARUDA_UTILS_CALLBACKS_HPP