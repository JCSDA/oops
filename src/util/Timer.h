/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef UTIL_TIMER_H_
#define UTIL_TIMER_H_

// #include <chrono>
#include <sys/time.h>

#include <string>

#include <boost/noncopyable.hpp>

namespace util {

// -----------------------------------------------------------------------------

class Timer : private boost::noncopyable {
 public:
  Timer(const std::string &, const std::string &);
  ~Timer();

 private:
//  const std::chrono::time_point<std::chrono::system_clock> start_;
  struct timeval start_;
  const std::string class_;
  const std::string method_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // UTIL_TIMER_H_
