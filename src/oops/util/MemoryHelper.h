/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <map>
#include <string>

#include "oops/util/Printable.h"

namespace util {

// -----------------------------------------------------------------------------

class MemoryHelper : public util::Printable {
 public:
  static void start();
  static void stop();
  static void add(const std::string &, const double);
  MemoryHelper(const MemoryHelper&) = delete;
  ~MemoryHelper();

 private:
  static MemoryHelper & getHelper();
  MemoryHelper();
  void print(std::ostream &) const;

  bool on_;
  std::map<std::string, std::array<double, 3>> stats_;
};

// -----------------------------------------------------------------------------

}  // namespace util

