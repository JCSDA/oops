/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_TIMERHELPER_H_
#define OOPS_UTIL_TIMERHELPER_H_

#include <iostream>
#include <map>
#include <memory>
#include <string>
// #include <chrono>

#include <boost/noncopyable.hpp>
#include "eckit/mpi/Comm.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Printable.h"
#include "oops/util/TimerTree.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class Timer;

// -----------------------------------------------------------------------------

class TimerHelper : public util::Printable,
                    private boost::noncopyable {
 public:
  static void start(const std::string &, const std::string &, const eckit::Configuration &);
  static void stop();
  static void init(const std::string &);
  static void add(const std::string &, const double, const bool);
//               const std::chrono::duration<double> &);
  static void setComm(const eckit::mpi::Comm & comm);
  ~TimerHelper();

 private:
  static TimerHelper & getHelper();
  TimerHelper();
  void print(std::ostream &) const;

  bool on_;
  std::map< std::string, double > timers_;
  std::map< std::string, int > counts_;
  std::unique_ptr<Timer> total_;
  std::string totalname_;
  const eckit::mpi::Comm * comm_;
  bool hierarchical_;
  std::unique_ptr<TimerTree> root_;
  TimerTree * current_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_TIMERHELPER_H_
