/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef UTIL_TIMERHELPER_H_
#define UTIL_TIMERHELPER_H_

// #include <chrono>
#include <iostream>
#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include "util/Printable.h"

namespace util {
  class Timer;

// -----------------------------------------------------------------------------

class TimerHelper : public util::Printable,
                    private boost::noncopyable {
 public:
  static void start();
  static void stop();
  static void add(const std::string &, const double);
//               const std::chrono::duration<double> &);
  ~TimerHelper();

 private:
  static TimerHelper & getHelper();
  TimerHelper();
  void print(std::ostream &) const;

  bool on_;
  std::map< std::string, double > timers_;
  std::map< std::string, int > counts_;
  boost::scoped_ptr<Timer> total_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // UTIL_TIMERHELPER_H_
