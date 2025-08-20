/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/gatherPrint.h"

namespace util {

// -----------------------------------------------------------------------------
void gatherAndPrint(std::ostream & os, const std::vector<char> & buffer,
                    const eckit::mpi::Comm & comm, size_t maxlen) {
  if (comm.size() > 1) {
    const size_t fulllen = buffer.size();
    ASSERT(fulllen % maxlen == 0);
    const size_t num_parts = fulllen / maxlen;
    std::vector<char> vglob(fulllen * comm.size());

    comm.gather(buffer, vglob, 0);

    if (comm.rank() == 0) {
      auto it = vglob.begin();
      for (size_t jpe = 0; jpe < comm.size(); ++jpe) {
        for (size_t jpart = 0; jpart < num_parts; ++jpart) {
          // Extract part of the buffer
          std::string spe(it, it + maxlen);
          std::size_t last = spe.find_last_not_of('#');
          if (last != std::string::npos) spe.erase(last + 1);
          os << spe;
          it += maxlen;
        }
      }
    }
  } else {
    // For single process, just convert buffer to string and print
    std::string str(buffer.begin(), buffer.end());
    std::size_t last = str.find_last_not_of('#');
    if (last != std::string::npos) str.erase(last+1);
    os << str;
  }
}

}  // namespace util
