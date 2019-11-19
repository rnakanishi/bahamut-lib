#include <utils/statistics.hpp>

#include <fstream>
#include <sstream>
#include <iostream>

namespace Ramuh {

Statistics::Statistics() {}

void Statistics::addDoubleField(std::string field) {
  _doubleComponents[field] = std::vector<double>();
}

void Statistics::addIntField(std::string field) {
  _intComponents[field] = std::vector<int>();
}

void Statistics::registerComponent(std::string field, double value) {
  auto it = _doubleComponents.find(field);
  if (it == _doubleComponents.end()) {
    _doubleComponents[field] = std::vector<double>();
  }
  _doubleComponents[field].emplace_back(value);
}

void Statistics::registerComponent(std::string field, int value) {
  auto it = _intComponents.find(field);
  if (it == _intComponents.end()) {
    _intComponents[field] = std::vector<int>();
  }
  _intComponents[field].emplace_back(value);
}

void Statistics::writeToFile(std::string filename) {
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

  // Components names: double fields and then int fields
  for (auto field : _doubleComponents) {
    buffer << field.first << ", ";
    steps = (steps < field.second.size()) ? field.second.size() : steps;
  }
  for (auto field : _intComponents) {
    buffer << field.first << ", ";
    steps = (steps < field.second.size()) ? field.second.size() : steps;
  }
  buffer.seekp(-2, buffer.cur);
  buffer << std::endl;
  file << buffer.str();

  // Corresponding values for double comps and int comps
  for (size_t i = 0; i < steps; i++) {
    buffer = std::ostringstream("");
    for (auto field : _doubleComponents) {
      if (field.second.size() >= steps - 1) {
        double timestamp;
        if (field.second.empty())
          timestamp = 0.0;
        else
          timestamp = field.second[i];
        buffer << timestamp << ", ";
      }
    }
    for (auto field : _intComponents) {
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
