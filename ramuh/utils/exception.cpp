#include <utils/exception.hpp>

namespace Ramuh {
Exception::Exception(const std::string &msg) : _message(msg) {}

const char *Exception::what() const noexcept { return _message.c_str(); }

} // namespace Ramuh
