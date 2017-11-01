/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "util/Timer.h"

// #include <chrono>
#include <string>
#include <sys/time.h>

#include <boost/noncopyable.hpp>

#include "util/TimerHelper.h"

namespace util {

// -----------------------------------------------------------------------------

Timer::Timer(const std::string & cl, const std::string & met)
    : class_(cl), method_(met) {
//      start_(std::chrono::high_resolution_clock::now()) {
  gettimeofday(&start_, NULL);
}

// -----------------------------------------------------------------------------

Timer::~Timer() {
//  end = std::chrono::high_resolution_clock::now();
//  std::chrono::duration<double> dt(end-start_);
//  timers_.addTime(class_, method_, dt);
  struct timeval end;
  gettimeofday(&end, NULL);

// get difference in milliseconds
  double dt = (end.tv_sec-start_.tv_sec)*1000.0
            + (end.tv_usec-start_.tv_usec)/1000.0;
  std::string name = class_ + "::" + method_;

  TimerHelper::add(name, dt);
};

// -----------------------------------------------------------------------------

}  // namespace util
