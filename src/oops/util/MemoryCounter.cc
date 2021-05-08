/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/MemoryCounter.h"

#include <string>

#include "eckit/system/ResourceUsage.h"
#include "oops/util/MemoryHelper.h"

namespace util {

// -----------------------------------------------------------------------------

MemoryCounter::MemoryCounter(const std::string & classname)
  : name_(classname), rss_(eckit::system::ResourceUsage().maxResidentSetSize())
{}

// -----------------------------------------------------------------------------

MemoryCounter::~MemoryCounter() {
  size_t current = eckit::system::ResourceUsage().maxResidentSetSize();
// convert bytes to MiB
  double mem = (static_cast<double>(current)-static_cast<double>(rss_))/1048576.0;
  MemoryHelper::add(name_, mem);
}

// -----------------------------------------------------------------------------

}  // namespace util

