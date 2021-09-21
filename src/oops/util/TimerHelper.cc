/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/TimerHelper.h"

#include <cmath>
#include <iomanip>
#include <string>

#include "eckit/io/Buffer.h"
#include "eckit/mpi/Comm.h"
#include "eckit/serialisation/ResizableMemoryStream.h"

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

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
  oops::Log::stats() << getHelper() << std::endl;
  getHelper().timers_.clear();
  getHelper().counts_.clear();
  getHelper().on_ = false;
}

// -----------------------------------------------------------------------------

void TimerHelper::add(const std::string & name, const double dt) {
  if (getHelper().on_) {
    getHelper().timers_[name] += dt;
    getHelper().counts_[name] += 1;
  }
}

// -----------------------------------------------------------------------------

TimerHelper::TimerHelper(): on_(false), timers_(), counts_(), total_() {}

// -----------------------------------------------------------------------------

TimerHelper::~TimerHelper() {}

// -----------------------------------------------------------------------------

void TimerHelper::print(std::ostream & os) const {
  typedef std::map<std::string, double>::const_iterator cit;

  {
    // Local timing statistics
    int table_width = 92;
    os << " " << std::endl;
    os << std::string(table_width, '-') << std::endl;
    std::string title = " Timing Statistics ";
    float title_half_width = (table_width-title.size())/2.;
    os << std::string(std::floor(title_half_width), '-')
       << title << std::string(std::ceil(title_half_width), '-') << std::endl
       << std::string(table_width, '-') << std::endl
       << std::setw(52) << std::left << "Name " << ": "
       << std::setw(12) << std::right << "total (ms)"
       << std::setw(8) << std::right << "count"
       << std::setw(18) << std::right << "time/call (ms)" << std::endl;
    for (cit jt = timers_.begin(); jt != timers_.end(); ++jt) {
      int icount = counts_.at(jt->first);
      os << std::setw(52) << std::left << jt->first
         << ": " << std::setw(12) << std::right << std::fixed << std::setprecision(2) << jt->second
         << std::setw(8) << icount
         << std::setw(18) << std::right << std::fixed << std::setprecision(4) << jt->second/icount
         << std::endl;
    }
    os << std::string(std::floor(title_half_width), '-')
       << title << std::string(std::ceil(title_half_width), '-') << std::endl;
  }
// For MPI applications, gather and print statistics across tasks
  size_t ntasks = oops::mpi::world().size();
  if (ntasks > 1) {
    int tag = 1234;
//  Tasks send their numbers to task 0
    if (oops::mpi::world().rank() > 0) {
      eckit::Buffer bufr(8000);
      eckit::ResizableMemoryStream sstr(bufr);
      sstr << timers_.size();
      for (cit jt = timers_.begin(); jt != timers_.end(); ++jt) {
        sstr << jt->first << jt->second;
      }
      oops::mpi::world().send(static_cast<const char*>(bufr.data()), sstr.position(), 0, tag);
    } else {  // Task 0
//    Structure for global statistics
      std::map<std::string, std::array<double, 3>> stats;
      for (cit jt = timers_.begin(); jt != timers_.end(); ++jt) {
        stats[jt->first].fill(jt->second);
      }
//    Task 0 receives stats from other tasks
      for (size_t from = 1; from < ntasks; ++from) {
        eckit::mpi::Status st = oops::mpi::world().probe(from, tag);
        size_t size = oops::mpi::world().getCount<char>(st);
        eckit::Buffer bufr(size);
        bufr.zero();

        oops::mpi::world().receive(static_cast<char*>(bufr.data()), bufr.size(), from, tag);
        eckit::ResizableMemoryStream sstr(bufr);

        std::string name;
        double time;
        size_t ntimes;
        sstr >> ntimes;
        for (size_t jj = 0; jj < ntimes; ++jj) {
          sstr >> name;
          sstr >> time;
          std::map<std::string, std::array<double, 3>>::iterator it = stats.find(name);
          if (it == stats.end()) {
            stats[name].fill(time);
          } else {
            if (time < it->second[0]) it->second[0] = time;
            if (time > it->second[1]) it->second[1] = time;
            it->second[2] += time;
          }
        }
      }
//    Print global statistics
      int table_width = 114;
      std::ostringstream title_s;
      title_s << " Parallel Timing Statistics (" << std::setw(4) << ntasks << " MPI tasks) ";
      std::string title = title_s.str();
      float title_half_width = (table_width-title.size())/2.;
      os << std::endl << std::string(table_width, '-') << std::endl
         << std::string(std::floor(title_half_width), '-')
         << title << std::string(std::ceil(title_half_width), '-') << std::endl
         << std::string(table_width, '-') << std::endl
         << std::setw(52) << std::left << "Name " << ": "
         << std::setw(12) << std::right << "min (ms)"
         << std::setw(12) << std::right << "max (ms)"
         << std::setw(12) << std::right << "avg (ms)"
         << std::setw(12) << std::right << "% total"
         << std::setw(12) << std::right << "imbal (%)"
         << std::endl;
      double total = stats["util::Timers::Total"][2]/ntasks;
      stats["util::Timers::measured"].fill(0.0);
      for (std::map<std::string, std::array<double, 3>>::iterator jt = stats.begin();
           jt != stats.end(); ++jt) {
//      Count what has been measured (to be compared with total run time)
        if (jt->first.substr(0, 4) == "oops") {
          for (size_t jj = 0; jj < 3; ++jj) {
            stats["util::Timers::measured"][jj] += jt->second[jj];
          }
        }
//      Only print for contributions greater than 0.1% of total
        double avg = jt->second[2]/ntasks;
        if (avg / total > 0.001)
          os << std::setw(52) << std::left << jt->first << ": "
             << std::setw(12) << std::right << std::fixed << std::setprecision(2)
             << std::setw(12) << jt->second[0]
             << std::setw(12) << jt->second[1]
             << std::setw(12) << avg
             << std::setw(12) << avg / total * 100.0
             << std::setw(12) << (jt->second[1] - jt->second[0]) / avg * 100.0
             << std::endl;
      }
      os << std::string(std::floor(title_half_width), '-')
         << title << std::string(std::ceil(title_half_width), '-') << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
