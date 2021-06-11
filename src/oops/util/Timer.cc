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

#include <iomanip>
#include <ratio>

#include "oops/util/Logger.h"
#include "oops/util/TimerHelper.h"

namespace util {

namespace {  // Local global constants
  // Width to print variable names
  constexpr int METHOD_PRINT_WIDTH = 60;
}

// -----------------------------------------------------------------------------

/** InitTime
 *
 * A class to be static initialized that measures time deltas since program initialization
 */
class InitTime
{
 public:
  InitTime() : init_t_(Timer::ClockT::now())
  { }

  double elapsed_sec(const Timer::TimeT &t) {
    return std::chrono::duration<double>(t - init_t_).count();
  }

 private:
  Timer::TimeT init_t_;
};

static InitTime init_time;  // static instance representing program static init time

// -----------------------------------------------------------------------------

Timer::Timer(const std::string & class_name, const std::string & method_name)
    : name_(class_name + "::" + method_name), start_(ClockT::now())
{ }

// -----------------------------------------------------------------------------

Timer::~Timer() {
  std::chrono::duration<double, std::milli> dt = ClockT::now() - start_;  // elapsed millisecs
  TimerHelper::add(name_, dt.count());
}

// -----------------------------------------------------------------------------

LoggingTimer::LoggingTimer(const std::string & class_name, const std::string & method_name)
    : LoggingTimer(class_name, method_name, oops::Log::timer())
{ }

// -----------------------------------------------------------------------------

LoggingTimer::LoggingTimer(const std::string & class_name,
                           const std::string & method_name, std::ostream& log)
    : Timer(class_name, method_name), log_(log)
{
  double st = init_time.elapsed_sec(start_);
  log_ << std::left << std::setw(METHOD_PRINT_WIDTH) << std::setfill('.') << name_
       << " Start: " << std::fixed << std::setprecision(2) << st << std::endl;
}

// -----------------------------------------------------------------------------

LoggingTimer::~LoggingTimer()
{
  double st = init_time.elapsed_sec(start_);
  double et = init_time.elapsed_sec(ClockT::now());
  log_ << std::left << std::setw(METHOD_PRINT_WIDTH) << std::setfill('.') << name_
       << ".. End: " << std::fixed << std::setprecision(2) << et
       << " Elapsed: " << et-st << " sec" << std::endl;
}

// -----------------------------------------------------------------------------

double timeStamp() {
  return init_time.elapsed_sec(Timer::ClockT::now());
}

// -----------------------------------------------------------------------------

}  // namespace util
