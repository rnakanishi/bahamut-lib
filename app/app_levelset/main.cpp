#include "levelset_fluid.h"
#include <iostream>

int main(void) {
  LevelSetFluid sim;

  std::cerr << "Initialized grid with size " << sim.size()
            << "  and resolution " << sim.h() << std::endl;

  return 0;
}