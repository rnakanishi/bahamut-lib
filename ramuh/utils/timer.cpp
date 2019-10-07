#include <utils/timer.hpp>
#include <iostream>
#include <iomanip>

namespace Ramuh {

Timer::Timer() {
  _start = _lastLap = _end = std::chrono::steady_clock::now();
  _longestName = 0;
  std::cout << "\033[21;32m===[Timer created]=== \033[0m\n";
}

void Timer::reset() { _start = _lastLap = std::chrono::steady_clock::now(); }

void Timer::clearAll() {
  for (auto &comp : _components) {
    comp.second = 0.0;
  }
}

void Timer::lap() { _lastLap = std::chrono::steady_clock::now(); }

double Timer::getEllapsedTime() {
  auto timeNow = std::chrono::steady_clock::now();
  std::chrono::duration<double> duration = timeNow - _lastLap;
  _lastLap = std::chrono::steady_clock::now();
  return duration.count();
}

double Timer::registerTime(const std::string &name, bool silence) {
  auto it = _components.find(name);
  double timeInterval = getEllapsedTime();
  if (!silence)
    std::cout << "\033[21;32m[TIMER]: \033[0m" << name << " took "
              << timeInterval << std::endl;
  if (it == _components.end())
    _components[name] = 0.0;
  _components[name] += timeInterval;
  if (name.length() > _longestName)
    _longestName = name.length();
  return 0.0; // FIXME: registeTime: fix return value
}

double Timer::getTotalTime() {
  auto timeNow = std::chrono::steady_clock::now();
  std::chrono::duration<double> duration = timeNow - _start;
  return duration.count();
}

void Timer::evaluateComponentsTime() {
  double total = getTotalTime();
  std::cout << "===============================\n";
  for (auto comp : _components) {
    std::cout << std::setw(_longestName) << comp.first << ":\t"
              << std::setprecision(4) << comp.second << "\t(" << std::setw(4)
              << std::setprecision(4) << comp.second / total * 100 << "\%)\n";
  }
  std::cout << "Total time: " << getTotalTime() << std::endl;
  std::cout << "===============================\n";
}

} // namespace Ramuh