/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Baudouin Raoult
/// @author Tiago Quintino
/// @date   December 2016

#include <algorithm>
#include <string>

#include "eckit/log/OStreamTarget.h"
#include "eckit/log/PrefixTarget.h"
#include "eckit/utils/Translator.h"

#include "oops/util/LibOOPS.h"

namespace oops {

//------------------------------------------------------------------------------

static LibOOPS liboops;

LibOOPS::LibOOPS() : Library("oops") {
  const char* e = ::getenv("OOPS_TRACE");
  if (e) {
    trace_ = eckit::Translator<std::string, bool>()(e);
  }
}

LibOOPS::~LibOOPS() {
}

LibOOPS& LibOOPS::instance() {
  return liboops;
}

void LibOOPS::finalise() {
  eckit::Log::flush();

  // Make sure that these specialised channels that wrap eckit::Log::info() are
  // destroyed before eckit::Log::info gets destroyed.
  // Just in case someone still tries to log, we reset to empty channels.
  traceChannel_.reset(new eckit::Channel());
  statsChannel_.reset(new eckit::Channel());
  testChannel_. reset(new eckit::Channel());
}

const void* LibOOPS::addr() const { return this; }

std::string LibOOPS::version() const { return "0.0.0"; }

std::string LibOOPS::gitsha1(unsigned int count) const {
  std::string sha1("");
  if (sha1.empty()) {
    return "not available";
  }
  return sha1.substr(0, std::min(count, 40u));
}

eckit::Channel& LibOOPS::traceChannel() const {
  if (traceChannel_) { return *traceChannel_; }
  if (trace_) {
    traceChannel_.reset(new eckit::Channel(
      new eckit::PrefixTarget("OOPS_TRACE", new eckit::OStreamTarget(eckit::Log::info()))));
  } else {
    traceChannel_.reset(new eckit::Channel());
  }
  return *traceChannel_;
}

eckit::Channel& LibOOPS::statsChannel() const {
  if (statsChannel_) { return *statsChannel_; }
  statsChannel_.reset(new eckit::Channel(
    new eckit::PrefixTarget("OOPS_STATS", new eckit::OStreamTarget(eckit::Log::info()))));
  return *statsChannel_;
}

eckit::Channel& LibOOPS::testChannel() const {
  if (testChannel_) { return *testChannel_; }
  testChannel_.reset(new eckit::Channel(
    new eckit::PrefixTarget("Test     :", new eckit::OStreamTarget(eckit::Log::info()))));
  return *testChannel_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

