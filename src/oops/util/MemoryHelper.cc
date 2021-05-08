/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/MemoryHelper.h"

#include <cmath>
#include <iomanip>
#include <string>

#include "eckit/io/ResizableBuffer.h"
#include "eckit/mpi/Comm.h"
#include "eckit/serialisation/ResizableMemoryStream.h"
#include "eckit/system/ResourceUsage.h"

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/MemoryCounter.h"

namespace util {

// -----------------------------------------------------------------------------

MemoryHelper & MemoryHelper::getHelper() {
  static MemoryHelper theHelper;
  return theHelper;
}

// -----------------------------------------------------------------------------

void MemoryHelper::start() {
  getHelper().on_ = true;
}

// -----------------------------------------------------------------------------

void MemoryHelper::stop() {
  oops::Log::stats() << getHelper() << std::endl;
  getHelper().stats_.clear();
  getHelper().on_ = false;
}

// -----------------------------------------------------------------------------

void MemoryHelper::add(const std::string & name, const double rss) {
  if (getHelper().on_) {
    getHelper().stats_[name][0] += 1.0;
    getHelper().stats_[name][1] += rss;
    if (rss > getHelper().stats_[name][2]) getHelper().stats_[name][2] = rss;
  }
}

// -----------------------------------------------------------------------------

MemoryHelper::MemoryHelper(): on_(false), stats_() {}

// -----------------------------------------------------------------------------

MemoryHelper::~MemoryHelper() {}

// -----------------------------------------------------------------------------

void MemoryHelper::print(std::ostream & os) const {
  size_t ntasks = oops::mpi::world().size();
  int tag = 1235;
  size_t rssbyte = eckit::system::ResourceUsage().maxResidentSetSize();
  double rss = static_cast<double>(rssbyte)/(1024*1024);

  if (oops::mpi::world().rank() > 0) {
//  Tasks send their numbers to task 0
    eckit::ResizableBuffer bufr(8000);
    eckit::ResizableMemoryStream sstr(bufr);
    sstr << stats_.size();
    for (auto jt = stats_.begin(); jt != stats_.end(); ++jt) {
      sstr << jt->first << jt->second[0] << jt->second[1] << jt->second[2];
    }
    sstr << rss;
    oops::mpi::world().send(static_cast<const char*>(bufr.data()), sstr.position(), 0, tag);
  } else {
//  Task 0 receives stats from other tasks
    std::map<std::string, std::array<double, 3>> stats = stats_;
    double rssmin = rss;
    double rssmax = rss;
    double rsstot = rss;
    for (size_t from = 1; from < ntasks; ++from) {
      eckit::mpi::Status st = oops::mpi::world().probe(from, tag);
      size_t size = oops::mpi::world().getCount<char>(st);
      eckit::ResizableBuffer bufr(size);
      bufr.zero();

      oops::mpi::world().receive(static_cast<char*>(bufr.data()), bufr.size(), from, tag);
      eckit::ResizableMemoryStream sstr(bufr);

      std::string name;
      double count;
      double sum;
      double max;
      size_t nstats;
      sstr >> nstats;
      for (size_t jj = 0; jj < nstats; ++jj) {
        sstr >> name;
        sstr >> count;
        sstr >> sum;
        sstr >> max;
        stats[name][0] += count;
        stats[name][1] += sum;
        if (max > stats[name][2]) stats[name][2] = max;
      }
      sstr >> rss;
      if (rss < rssmin) rssmin = rss;
      if (rss > rssmax) rssmax = rss;
      rsstot += rss;
    }

//  Print global statistics
    int table_width = 64;
    std::ostringstream title_s;
    title_s << " Memory Statistics (" << std::setw(5) << ntasks << " MPI tasks) ";
    std::string title = title_s.str();
    float title_half_width = (table_width-title.size())/2.;
    os << std::endl << std::string(table_width, '-') << std::endl
       << std::string(std::floor(title_half_width), '-')
       << title << std::string(std::ceil(title_half_width), '-') << std::endl
       << std::string(table_width, '-') << std::endl
       << std::setw(30) << std::left << " "
       << std::setw(11) << std::right << "Count/task"
       << std::setw(11) << std::right << "Avg (MiB)"
       << std::setw(12) << std::right << "Max (MiB)"
       << std::endl;
    for (auto jt = stats.begin(); jt != stats.end(); ++jt) {
      os << std::setw(30) << std::left << jt->first << ": "
         << std::setw(12) << std::right << std::fixed << std::setprecision(0)
         << std::setw(8)  << jt->second[0]/ntasks
         << std::setw(12) << std::right << std::fixed << std::setprecision(2)
         << std::setw(12) << jt->second[1]/jt->second[0]
         << std::setw(12) << jt->second[2]
         << std::endl;
    }
    os << std::string(table_width, '-') << std::endl;

    os << "Total memory used per MPI task: min = "  << std::setprecision(0)
       << rssmin << " MiB, max = " << rssmax << " MiB." << std::endl;
    os << "Total memory used: "  << std::setprecision(0) << rsstot << " MiB." << std::endl;
    os << std::string(table_width, '-') << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
