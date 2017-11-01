/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "util/TimerHelper.h"

// #include <chrono>
#include <iomanip>
#include <string>

#include "util/Logger.h"
#include "util/Timer.h"

using oops::Log;

namespace util {

// -----------------------------------------------------------------------------

TimerHelper & TimerHelper::getHelper() {
  static TimerHelper theHelper;
  return theHelper;
}

// -----------------------------------------------------------------------------

void TimerHelper::start() {
  getHelper().on_ = true;
  getHelper().total_.reset(new Timer("util::Timers", "Total"));
}

// -----------------------------------------------------------------------------

void TimerHelper::stop() {
  getHelper().total_.reset();
  Log::stats() << getHelper() << std::endl;
  getHelper().timers_.clear();
  getHelper().counts_.clear();
  getHelper().on_ = false;
}

// -----------------------------------------------------------------------------

void TimerHelper::add(const std::string & name, const double dt) {
//                      const std::chrono::duration<double> & dt) {
  if (getHelper().on_) {
    getHelper().timers_[name] += dt;
    getHelper().counts_[name] += 1;
//    double secs = dt.count();
//    getHelper().timers_[name] += secs;
  }
}

// -----------------------------------------------------------------------------

TimerHelper::TimerHelper(): on_(false), timers_(), counts_(), total_() {}

// -----------------------------------------------------------------------------

TimerHelper::~TimerHelper() {}

// -----------------------------------------------------------------------------

void TimerHelper::print(std::ostream & os) const {
  os << " " << std::endl;
  os << "---------------------------------------------------------------------" << std::endl;
  os << "------------------------- Timing Statistics -------------------------" << std::endl;
  os << "---------------------------------------------------------------------" << std::endl;
  typedef std::map<std::string, double>::const_iterator cit;
  for (cit jt = timers_.begin(); jt != timers_.end(); ++jt) {
    int icount = counts_.at(jt->first);
    os << std::setw(52) << std::left << jt->first
       << ": " << std::setw(12) << std::right << jt->second << " ms"
       << "  " << std::setw(6) << icount << "   " << std::setw(12) << std::right << jt->second/icount << " ms/call"
       << std::endl;
  }
  os << "------------------------- Timing Statistics -------------------------" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace util
