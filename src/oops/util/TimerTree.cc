/*
 * (C) Copyright 2025-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
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

TimerTree * TimerTree::goUp() {
  return parent_;
}

// -----------------------------------------------------------------------------

void TimerTree::addTime(const double time) {
  time_+= time;
  ++calls_;
}

// -----------------------------------------------------------------------------

void TimerTree::printTree(std::ostream & os, const std::string & prefix, bool isLast) const {
  // Print current node
  os << prefix;
  os << (isLast ? "└── " : "├── ");
  os << name_ << " (" << std::fixed << std::setprecision(2) << time_<< " ms, "
     << calls_ << " calls)" << std::endl;

  // Print children
  if (!children_.empty()) {
    auto it = children_.begin();
    for (size_t i = 0; i < children_.size(); ++i, ++it) {
      bool childIsLast = (i == children_.size() - 1);
      std::string childPrefix = prefix + (isLast ? "    " : "│   ");
      it->second->printTree(os, childPrefix, childIsLast);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
