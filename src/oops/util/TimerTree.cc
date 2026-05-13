/*
 * (C) Copyright 2025-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "oops/util/TimerTree.h"

namespace util {

// -----------------------------------------------------------------------------

TimerTree::TimerTree(const std::string & name)
  : name_(name), time_(0.0), calls_(0), parent_(nullptr)
{}

// -----------------------------------------------------------------------------

TimerTree * TimerTree::goDown(const std::string & name) {
  if (children_.find(name) == children_.end()) {
    children_[name] = std::make_unique<TimerTree>(name);
    children_[name]->parent_ = this;
  }
  return children_[name].get();
}

// -----------------------------------------------------------------------------

TimerTree * TimerTree::goUp() const {
  return parent_;
}

// -----------------------------------------------------------------------------

void TimerTree::addTime(const double time) {
  time_+= time;
  ++calls_;
}

// -----------------------------------------------------------------------------

double TimerTree::time() const {
  return time_;
}

// -----------------------------------------------------------------------------

void TimerTree::printTree(std::ostream & os, const std::string & prefix, const bool isLast) const {
  // Print current node
  os << prefix;
  os << (isLast ? "└── " : "├── ");
  os << name_ << " (" << std::fixed << std::setprecision(2) << time_<< " ms, "
     << calls_ << " calls)" << std::endl;

  // Print children sorted in descending order of time taken.
  if (!children_.empty()) {
    std::vector<std::pair<std::string, double>> childTimes;
    childTimes.reserve(children_.size());
    for (auto it = children_.begin(); it != children_.end(); ++it) {
      childTimes.emplace_back(it->first, it->second->time());
    }
    std::sort(childTimes.begin(), childTimes.end(),
              [](auto &a, auto &b) {return a.second > b.second;});
    for (size_t i = 0; i < childTimes.size(); ++i) {
      const bool childIsLast = (i == childTimes.size() - 1);
      const std::string childPrefix = prefix + (isLast ? "    " : "│   ");
      children_.at(childTimes[i].first)->printTree(os, childPrefix, childIsLast);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
