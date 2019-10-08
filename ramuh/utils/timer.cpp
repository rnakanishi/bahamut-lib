#include <utils/timer.hpp>
#include <iostream>
#include <iomanip>

namespace Ramuh {

Timer::Timer() {
  _creation = _start = _lastLap = _end = std::chrono::steady_clock::now();
  _longestName = 0;
  _resetTimes = 1;
  _silence = false;
  std::cout << "\033[21;32m===[Timer created]=== \033[0m\n";
}

void Timer::setLogging(bool silence) { _silence = silence; }

void Timer::reset() {
  _start = _lastLap = std::chrono::steady_clock::now();
  _resetTimes++;
}

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
  if (it == _components.end()) {
    _components[name] = 0.0;
    _cumulative[name] = 0.0;
    _calls[name] = 0;
  }
  _components[name] += timeInterval;
  _cumulative[name] += timeInterval;
  _calls[name]++;
  if (name.length() > _longestName)
    _longestName = name.length();
  return 0.0; // FIXME: registeTime: fix return value
}

double Timer::getTotalTime() {
  auto timeNow = std::chrono::steady_clock::now();
  std::chrono::duration<double> duration = timeNow - _creation;
  return duration.count();
}

double Timer::getLapTime() {
  auto timeNow = std::chrono::steady_clock::now();
  std::chrono::duration<double> duration = timeNow - _start;
  return duration.count();
}

void Timer::evaluateComponentsTime() {
  double total = getLapTime();
  std::cout << "===============================\n";
  for (auto comp : _components) {
    std::cout << std::setw(_longestName) << comp.first << ":\t"
              << std::setprecision(4) << comp.second << "\t(" << std::setw(4)
              << std::setprecision(4) << comp.second / total * 100 << "\%)\n";
  }
  std::cout << "Total time: " << getLapTime() << std::endl;
  std::cout << "===============================\n";
}

void Timer::evaluateComponentsAverageTime() {
  double total = getTotalTime();
  std::cout << "=========AVERAGE TIME==========\n";
  for (auto comp : _cumulative) {
    int call = _calls[comp.first];

    std::cout << std::setw(_longestName) << comp.first << ":\t"
              << std::setprecision(4) << comp.second / call << "\t("
              << std::setw(4) << std::setprecision(4)
              << comp.second / total * 100 << "\%)\n";
  }
  std::cout << "Average time: " << getTotalTime() / _resetTimes << std::endl;
  std::cout << "===============================\n";
}

} // namespace Ramuh