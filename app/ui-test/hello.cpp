#include <iostream>
#include <gui/gui.hpp>


int main(int argc, char const *argv[]) {
  Garuda::GUI interface;

  interface.createWindow();

  interface.run();

  return 0;
}