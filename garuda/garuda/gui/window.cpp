#include <gui/window.hpp>

namespace Garuda {
Window::Window() {}

Window::Window(Eigen::Array2i size) : nanogui::Screen(size, "Garuda Window") {}
} // namespace Garuda
