#ifndef __RAMUH_UTILS_LOG_HPP__
#define __RAMUH_UTILS_LOG_HPP__

#include <string>

namespace Ramuh {
class Log {
public:
  Log(/* args */);
  ~Log();

  enum class ColorCode {
    FG_RED = 31,
    FG_GREEN = 32,
    FG_YELLOW = 33,
    FG_BLUE = 34,
    FG_MAGENTA = 35,
    FG_CYAN = 36,
    FG_WHITE = 37,
    FG_DEFAULT = 39
  };

  static void raiseError(const std::string &message);

  static void raiseWarning(const std::string &message);

  static void log(const std::string &message);

private:
  /* data */
};

} // namespace Ramuh

#endif