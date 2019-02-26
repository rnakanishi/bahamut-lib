#include <iostream>

#define NOT_IMPLEMENTED()                                                      \
  {                                                                            \
    std::cerr << "Not yet implemented: " << __FILE__ << ", " << __LINE__       \
              << std::endl;                                                    \
  }