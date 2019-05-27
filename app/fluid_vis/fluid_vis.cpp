#include "fluid_gui.hpp"

int main(int argc, char const *argv[]) {
  FluidGUI interface;

  interface.createWindow();
  interface.prepareCallbacks();
  interface.loadScene();

  interface.run();

  return 0;
}