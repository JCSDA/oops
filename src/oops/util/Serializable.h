/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

namespace util {

// -----------------------------------------------------------------------------

/// Add interfaces for methods to serialize objects for sending via MPI.

class Serializable {
 public:
  Serializable() {}
  virtual ~Serializable() {}

/// Size of serialized object (mostly so c++ can allocate memory for fortran implementations)
  virtual size_t serialSize() const = 0;

/// Serialize (appends data to vector argument)
  virtual void serialize(std::vector<double> &) const = 0;

/// Deserialize (using data starting from index given by second arg, update index to end
///              of current object for next one)
  virtual void deserialize(const std::vector<double> &, size_t &) = 0;
};

// -----------------------------------------------------------------------------

}  // namespace util
