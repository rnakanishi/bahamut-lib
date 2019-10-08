#ifndef __RAMUH_UTILS_TIMER_HPP
#define __RAMUH_UTILS_TIMER_HPP
#include <ctime>
#include <map>
#include <string>
#include <chrono>

namespace Ramuh {
class Timer {
public:
  Timer();

  void setLogging(bool silence = true);

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
   * @param [bool:false] silence if set, the method does not display current
   * evaluation time, but it is computed internally
   * @return double
   */
  double registerTime(const std::string &name, bool silence = false);

  /**
   * @brief Get the Total Time  pasted since object instantiation, or since last
   * reset call was made
   *
   * @return double
   */
  double getTotalTime();
  double getLapTime();

  /**
   * @brief Prints every componet taken time.
   * This method DOES NOT reset components timer.
   *
   */
  void evaluateComponentsTime();

  /**
   * @brief Prints all components average time, based on how many times %reset()
   * method was called.
   *
   */
  void evaluateComponentsAverageTime();

protected:
  std::chrono::time_point<std::chrono::steady_clock,
                          std::chrono::duration<double>>
      _start, _end, _lastLap, _creation;
  int _longestName;
  std::map<std::string, double> _components;
  std::map<std::string, double> _cumulative;
  std::map<std::string, int> _calls;
  int _resetTimes;
  bool _silence;
};

} // namespace Ramuh

#endif
