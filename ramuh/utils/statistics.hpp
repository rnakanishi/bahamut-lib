#ifndef __RAMUH_STATISTICS_HPP__
#define __RAMUH_STATISTICS_HPP__

#include <string>
#include <utils/timer.hpp>
#include <memory>

namespace Ramuh {
class Statistics {

public:
  Statistics();

  void addDoubleField(std::string field);

  void addIntField(std::string field);

  void registerComponent(std::string field, double value);

  void registerComponent(std::string field, int value);

  void writeToFile(std::string filename);

private:
  std::map<std::string, std::vector<double>> _doubleComponents;
  std::map<std::string, std::vector<int>> _intComponents;
};
} // namespace Ramuh

#endif