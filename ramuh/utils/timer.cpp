#include <utils/timer.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

namespace Ramuh {

Timer::Timer() {
  _creation = _start = _lastLap = _end = std::chrono::steady_clock::now();
  _longestName = 0;
  _resetTimes = 1;
  _silence = false;
  std::cout << "\033[21;32m===[Timer created]=== \033[0m\n";

  _log["total"] = std::vector<double>();
}

void Timer::setLogging(bool silence) { _silence = silence; }

void Timer::reset() {
  _start = _lastLap = std::chrono::steady_clock::now();
  _resetTimes++;
  clearAll();
}

void Timer::clearAll() {
  for (auto &comp : _components) {
    comp.second = 0.0;
    _callsBeforeReset[comp.first] = 0;
  }
}

void Timer::lap() { _lastLap = std::chrono::steady_clock::now(); }

double Timer::getComponentTime(const std::string &name) {
  auto it = _components.find(name);
  if (it == _components.end()) {
    std::cout << "\033[21;33m[Warning]: \033[0mTimer['" << name
              << "'] not found.\n";
    return 0.0;
  }
  return it->second;
}

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
    _log[name] = std::vector<double>();
    _calls[name] = 0;
    _callsBeforeReset[name] = 0;
  }
  _components[name] += timeInterval;
  _cumulative[name] += timeInterval;
  _calls[name]++;
  _callsBeforeReset[name]++;
  if (name.length() > _longestName)
    _longestName = name.length();

  return timeInterval;
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

    int ncalls = _callsBeforeReset[comp.first];
    _log[comp.first].emplace_back(comp.second / (double)ncalls);
  }
  _log["total"].emplace_back(total);
  _components["total"] = total;
  std::cout << "Total time: " << total << std::endl;
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

void Timer::logToFile(std::string filename) {
  static int count = 0;
  std::ofstream file;

  file.open(filename, std::ofstream::out);
  if (!file.is_open()) {
    std::cerr << "\033[1;31mError\033[0m Timer::logToFile: Failed to open "
              << filename << std::endl;
    return;
  }

  std::ostringstream buffer;
  int steps = 0;
  for (auto field : _log) {
    buffer << field.first << ", ";
    steps = (steps < field.second.size()) ? field.second.size() : steps;
  }
  buffer.seekp(-2, buffer.cur);
  buffer << std::endl;
  file << buffer.str();

  for (size_t i = 0; i < steps; i++) {
    buffer = std::ostringstream("");
    for (auto field : _log) {
      if (field.second.size() >= steps - 1) {
        double timestamp;
        if (field.second.empty())
          timestamp = 0.0;
        else
          timestamp = field.second[i];
        buffer << timestamp << ", ";
      }
    }
    buffer.seekp(-2, buffer.cur);
    buffer << std::endl;
    file << buffer.str();
  }
  file.close();
  std::cout << "File written: " << filename << std::endl;
}

} // namespace Ramuh