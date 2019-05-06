#include <utils/timer.hpp>
#include <iostream>
#include <iomanip>

namespace Ramuh {

Timer::Timer() {
  _start = _lastLap = _end = std::clock();
  _longestName = 0;
}

void Timer::reset() { _start = _lastLap = std::clock(); }

void Timer::clearAll() {
  for (auto &comp : _components) {
    comp.second = 0.0;
  }
}

void Timer::lap() { _lastLap = std::clock(); }

double Timer::getEllapsedTime() {
  double seconds =
      (std::clock() - _lastLap) / static_cast<double>(CLOCKS_PER_SEC);
  _lastLap = std::clock();
  return seconds;
}

double Timer::registerTime(const std::string &name) {
  auto it = _components.find(name);
  if (it == _components.end())
    _components[name] = 0.0;
  _components[name] += getEllapsedTime();
  if (name.length() > _longestName)
    _longestName = name.length();
}

double Timer::getTotalTime() {
  return (std::clock() - _start) / static_cast<double>(CLOCKS_PER_SEC);
}

void Timer::evaluateComponentsTime() {
  double total = getTotalTime();
  std::cout << "===============================\n";
  for (auto comp : _components) {
    std::cout << std::setw(_longestName) << comp.first << ":\t"
              << std::setprecision(4) << comp.second << "\t(" << std::setw(4)
              << std::setprecision(4) << comp.second / total * 100 << "\%)\n";
  }
  std::cout << "Total time: " << getTotalTime << std::endl;
  std::cout << "===============================\n";
}

} // namespace Ramuh