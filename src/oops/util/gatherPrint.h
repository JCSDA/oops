/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"

namespace util {

// -----------------------------------------------------------------------------
/// Serializes an object into an existing string buffer with padding
/// for MPI communication
template <typename T>
void serializeForGather(std::vector<char> & buffer, const T & obj, std::ostream & os,
                        size_t maxlen = 10000) {
  std::stringstream ss;
  ss.setf(os.flags());
  ss.precision(os.precision());
  ss << obj;
  std::string sloc = ss.str();

  size_t start_pos = buffer.size();
  buffer.resize(buffer.size() + maxlen);

  // Copy string to buffer
  std::copy(sloc.begin(), sloc.end(), buffer.begin() + start_pos);

  // Pad remaining space with '#'
  std::fill(buffer.begin() + start_pos + sloc.size(), buffer.begin() + start_pos + maxlen, '#');
}

// -----------------------------------------------------------------------------

/// Gathers character buffers from all tasks and prints them in order on rank 0
void gatherAndPrint(std::ostream & os, const std::vector<char> & buffer,
                    const eckit::mpi::Comm & comm, size_t maxlen = 10000);

// -----------------------------------------------------------------------------

/// Collects prints from all the tasks in communicator to print in reproducible order
/// on task 0.
template <typename T>
void gatherPrint(std::ostream & os, const T & obj, const eckit::mpi::Comm & comm) {
  if (comm.size() > 1) {
    std::vector<char> buffer;
    serializeForGather(buffer, obj, os);
    gatherAndPrint(os, buffer, comm);
  } else {
    os << obj;
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
