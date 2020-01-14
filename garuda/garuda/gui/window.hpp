#ifndef __GARUDA_GUI_WINDOW_HPP__
#define __GARUDA_GUI_WINDOW_HPP__

#include <nanogui/nanogui.h>
#include <Eigen/Dense>

namespace Garuda {
class Window : public nanogui::Screen {
public:
  Window();

  Window(Eigen::Array2i size);

private:
};
}

#endif