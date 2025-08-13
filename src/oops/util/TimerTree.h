#ifndef OOPS_UTIL_TIMER_TREE_H_
#define OOPS_UTIL_TIMER_TREE_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

#include "oops/util/Timer.h"

namespace util {

class TimerTree {
 public:
  TimerTree(const std::string & name): name_(name), time_count_(0), called_count_(0), parent_(nullptr) {}
  ~TimerTree() {
    for (std::map<std::string, TimerTree*>::iterator it = children_.begin(); it != children_.end(); ++it) {
      delete it->second;
    }
  }
  TimerTree* goDown(const std::string & name) {
    if (children_.find(name) == children_.end()) {
      children_[name] = new TimerTree(name);
      children_[name]->parent_ = this;
    }
    return children_[name];
  }
  TimerTree* goUp() {
    return parent_;
  }
  void addTime(const int time) {
    time_count_ += time;
    ++called_count_;
  }
  
  // Print hierarchical tree structure like Unix 'tree' command
  void printTree(std::ostream& os, const std::string& prefix = "", bool isLast = true) const {
    // Print current node
    os << prefix;
    os << (isLast ? "└── " : "├── ");
    os << name_ << " (" << std::fixed << std::setprecision(2) << time_count_ << " ms, " 
       << called_count_ << " calls)" << std::endl;
    
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
  
 private:
  std::string name_;
  int time_count_, called_count_;
  TimerTree* parent_;
  std::map<std::string, TimerTree*> children_;
};
}





#endif  // OOPS_UTIL_TIMER_TREE_H_