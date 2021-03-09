/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_TIMER_H_
#define OOPS_UTIL_TIMER_H_

#include <chrono>
#include <ostream>
#include <string>

namespace util {

// -----------------------------------------------------------------------------

class Timer {
 public:
  using ClockT = std::chrono::steady_clock;  // Monotonic clock type
  using TimeT = ClockT::time_point;

  Timer(const std::string &class_name, const std::string &method_name);
  ~Timer();

  Timer(const Timer&) = delete;  // Non-copyable

 protected:
  std::string name_;
  TimeT start_;
};

class LoggingTimer : public Timer {
 public:
  LoggingTimer(const std::string &class_name, const std::string &method_name);
  LoggingTimer(const std::string &class_name, const std::string &method_name, std::ostream &log);
  ~LoggingTimer();

 protected:
  std::ostream& log_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_TIMER_H_
