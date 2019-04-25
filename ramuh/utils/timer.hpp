#ifndef __RAMUH_UTILS_TIMER_HPP
#define __RAMUH_UTILS_TIMER_HPP
#include <ctime>
#include <map>
#include <string>

namespace Ramuh {
class Timer {
public:
  Timer();

  void reset();

  void lap();

  double getEllapsedTime();

  double registerTime(const std::string &name);

  double getTotalTime();

  void evaluateComponentsTime();

protected:
  std::clock_t _start, _end, _lastLap;
  int _longestName;
  std::map<std::string, double> _components;
};

} // namespace Ramuh

#endif
