#ifndef __RAMUH_UTILS_EXCEPTION__
#define __RAMUH_UTILS_EXCEPTION__
#include <exception>
#include <string>

namespace Ramuh {

class Exception : public std::exception {
public:
  Exception(const std::string &msg);

  virtual const char *what() const noexcept override;

private:
  std::string _message;
};
} // namespace Ramuh

#endif
