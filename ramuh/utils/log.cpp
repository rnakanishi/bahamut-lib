#include <utils/log.hpp>
#include <iostream>

namespace Ramuh {

Log::Log(/* args */) {}

Log::~Log() {}

void Log::raiseError(const std::string &message) {
  std::cerr << "\033[1;31m[ERROR]\033[0m" << message << std::endl;
}

void Log::raiseWarning(const std::string &message) {
  std::cerr << "\033[21;33m[WARNING]\033[0m" << message << std::endl;
}

void Log::log(const std::string &message) {
  std::cerr << "\033[1;32m[LOG]\033[0m" << message << std::endl;
}

} // namespace Ramuh
