#ifndef __RAMUH_UTILS_TIMER_HPP
#define __RAMUH_UTILS_TIMER_HPP
#include <ctime>
#include <map>
#include <string>

namespace Ramuh {
class Timer {
public:
  Timer();

  /**
   * @brief Reset the timer and lap time, setting both to current time (zero)
   *
   */
  void reset();

  /**
   * @brief Clear all labels registered before the method call
   *
   */
  void clearAll();

  void lap();

  /**
   * @brief Get the time pasted since this method last call, ot last lap()
   * method call.
   *
   * @return double
   */
  double getEllapsedTime();

  /**
   * @brief Register a new label to evaluate. If the label already exists, then
   * the paste time is summed up.
   *  Prints time pasted since last evaluation, or
   * last label registered
   *
   * @param name
   * @return double
   */
  double registerTime(const std::string &name);

  /**
   * @brief Get the Total Time  pasted since object instantiation, or since last
   * reset call was made
   *
   * @return double
   */
  double getTotalTime();

  /**
   * @brief Prints every componet taken time.
   * This method DOES NOT reset components timer.
   *
   */
  void evaluateComponentsTime();

protected:
  std::clock_t _start, _end, _lastLap;
  int _longestName;
  std::map<std::string, double> _components;
};

} // namespace Ramuh

#endif
