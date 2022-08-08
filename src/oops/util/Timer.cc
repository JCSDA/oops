/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/Timer.h"

#include <chrono>

#include "oops/util/TimerHelper.h"

namespace util {

// -----------------------------------------------------------------------------

static std::chrono::steady_clock::time_point start_time(std::chrono::steady_clock::now());

static int nested_timers = 0;  // Only non-nested timers count towards measured time

// -----------------------------------------------------------------------------

Timer::Timer(const std::string & class_name, const std::string & method_name)
    : name_(class_name + "::" + method_name), start_()
{
  ++nested_timers;
  start_ = ClockT::now();
}

// -----------------------------------------------------------------------------

Timer::~Timer() {
  std::chrono::duration<double, std::milli> dt = ClockT::now() - start_;  // elapsed millisecs
  --nested_timers;
  // A top-level timer is created (when nested_timers == 0) in TimerHelper::start() for total time.
  // To count measured time (and establish timer coverage), we sum times from the timers 1 level
  // below this top-level timer. More-deeply nested timers would duplicate time if included.
  const bool include_timer_in_sum = (nested_timers == 1);
  TimerHelper::add(name_, dt.count(), include_timer_in_sum);
}

// -----------------------------------------------------------------------------

double timeStamp() {
  const std::chrono::steady_clock::time_point t(std::chrono::steady_clock::now());
  return std::chrono::duration<double>(t - start_time).count();
}

// -----------------------------------------------------------------------------

}  // namespace util
