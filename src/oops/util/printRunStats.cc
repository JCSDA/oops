/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/printRunStats.h"

#include <iomanip>
#include <string>
#include <vector>

#include "eckit/system/ResourceUsage.h"

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace util {

// -----------------------------------------------------------------------------

void printRunStats(const std::string & name, const bool alltasks,
                   const eckit::mpi::Comm & comm) {
  size_t rssbyte = eckit::system::ResourceUsage().maxResidentSetSize();
  double rss = static_cast<double>(rssbyte);

  double factor = 1.0e+6;
  std::string unit = " Mb";

  if (alltasks) {
    size_t ntasks = comm.size();
    std::vector<double> zss(ntasks);
    comm.gather(rss, zss, 0);

    if (comm.rank() == 0) {
      double rssmin = rss;
      double rssmax = rss;
      double rsstot = rss;
      for (size_t jj = 1; jj < ntasks; ++jj) {
        if (zss[jj] < rssmin) rssmin = zss[jj];
        if (zss[jj] > rssmax) rssmax = zss[jj];
        rsstot += zss[jj];
      }
      if (rsstot >= 1.0e+9) {factor = 1.0e+9; unit = " Gb";}
      oops::Log::stats() << std::left << std::setw(40) << name << " - Runtime: "
                         << std::fixed << std::right << std::setprecision(2)
                         << std::setw(9) << timeStamp() << " sec,  Memory: total: "
                         << std::setw(8) << rsstot / factor << unit;
      if (rssmin >= 1.0e+9) {factor = 1.0e+9; unit = " Gb";} else {factor = 1.0e+6; unit = " Mb";}
      oops::Log::stats() << ", per task: min = "
                         << std::right << std::setprecision(2) << std::fixed
                         << std::setw(8) << rssmin / factor << unit << ", max = "
                         << std::setw(8) << rssmax / factor << unit << std::endl;
    }
  } else {
    if (rss >= 1.0e+9) {factor = 1.0e+9; unit = " Gb";}
    oops::Log::stats() << std::left << std::setw(40) << name << " - Runtime: "
                       << std::fixed << std::right << std::setprecision(2)
                       << std::setw(8) << timeStamp() << " sec,  Local Memory: "
                       << std::setw(8) << rss / factor << unit << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
