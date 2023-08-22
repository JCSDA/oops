/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"

namespace util {

// -----------------------------------------------------------------------------

/// Collects prints from all the tasks in communicator to print in reproducible order
/// on task 0.
template <typename T>
void gatherPrint(std::ostream & os, const T & obj, const eckit::mpi::Comm & comm) {
  if (comm.size() > 1) {
    size_t maxlen = 10000;

    std::stringstream ss;
    ss.setf(os.flags());
    ss.precision(os.precision());
    ss << obj;
    std::string sloc = ss.str();
    std::vector<char> vloc(sloc.begin(), sloc.end());
    for (size_t jj = vloc.size(); jj < maxlen; ++jj) vloc.push_back('#');
    std::vector<char> vglob(maxlen*comm.size());

    comm.gather(vloc, vglob, 0);

    if (comm.rank() == 0) {
      auto it = vglob.begin();
      for (size_t jpe = 0; jpe < comm.size(); ++jpe) {
        std::string spe(it, it + maxlen);
        std::size_t last = spe.find_last_not_of('#');
        if (last != std::string::npos) spe.erase(last+1);
        os << spe;
        it += maxlen;
      }
    }
  } else {
    oops::Log::warning() << "No need to call gatherPrint";
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
