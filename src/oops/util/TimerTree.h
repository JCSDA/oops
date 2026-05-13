/*
 * (C) Copyright 2025-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace util {

class TimerTree {
 public:
  explicit TimerTree(const std::string &);
  ~TimerTree() = default;
  TimerTree * goDown(const std::string &);
  TimerTree * goUp() const;
  void addTime(const double);
  void printTree(std::ostream &, const std::string& prefix = "", const bool isLast = true) const;
  double time() const;

 private:
  std::string name_;
  double time_;
  int calls_;
  TimerTree * parent_;
  std::map<std::string, std::unique_ptr<TimerTree>> children_;
};

}  // namespace util

