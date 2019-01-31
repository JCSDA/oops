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

#include "oops/parallel/mpi/mpi.h"
#include "oops/util/LibOOPS.h"

namespace oops {

//------------------------------------------------------------------------------

bool getEnv(const std::string& env, bool default_value) {
  if (::getenv(env.c_str())) {return eckit::Translator<std::string, bool>()(::getenv(env.c_str()));}
  return default_value;
}

int getEnv(const std::string& env, int default_value) {
  if (::getenv(env.c_str())) {return eckit::Translator<std::string, int>()(::getenv(env.c_str()));}
  return default_value;
}

//------------------------------------------------------------------------------

static LibOOPS liboops;

LibOOPS::LibOOPS() : Library("oops"), rank_(oops::mpi::comm().rank()),
                     debug_(false), predebug_("OOPS_DEBUG"),
                     trace_(false), pretrace_("OOPS_TRACE") {
  const int it = getEnv("OOPS_TRACE", 0);
  if (it > 0 && rank_ == 0) trace_ = true;
  if (it < 0) {
    trace_ = true;
    pretrace_ += "[" + std::to_string(rank_) + "]";
  }
  const int id = getEnv("OOPS_DEBUG", 0);
  if (id > 0 && rank_ == 0) debug_ = true;
  if (id < 0) {
    debug_ = true;
    predebug_ += "[" + std::to_string(rank_) + "]";
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
  infoChannel_.reset(new eckit::Channel());
  debugChannel_.reset(new eckit::Channel());
  traceChannel_.reset(new eckit::Channel());
  statsChannel_.reset(new eckit::Channel());
  testChannel_. reset(new eckit::Channel());
}

const void* LibOOPS::addr() const {return this;}

std::string LibOOPS::version() const {return "0.0.0";}

std::string LibOOPS::gitsha1(unsigned int count) const {
  std::string sha1("");
  if (sha1.empty()) {
    return "not available";
  }
  return sha1.substr(0, std::min(count, 40u));
}

eckit::Channel& LibOOPS::traceChannel() const {
  if (traceChannel_) {return *traceChannel_;}
  if (trace_) {
    traceChannel_.reset(new eckit::Channel(
      new eckit::PrefixTarget(pretrace_, new eckit::OStreamTarget(eckit::Log::info()))));
  } else {
    traceChannel_.reset(new eckit::Channel());
  }
  return *traceChannel_;
}

eckit::Channel& LibOOPS::statsChannel() const {
  if (statsChannel_) {return *statsChannel_;}
  if (rank_ == 0) {
    statsChannel_.reset(new eckit::Channel(
      new eckit::PrefixTarget("OOPS_STATS", new eckit::OStreamTarget(eckit::Log::info()))));
  } else {
    statsChannel_.reset(new eckit::Channel());
  }
  return *statsChannel_;
}

eckit::Channel& LibOOPS::testChannel() const {
  if (testChannel_) {return *testChannel_;}
  if (rank_ == 0) {
    testChannel_.reset(new eckit::Channel(
      new eckit::PrefixTarget("Test     :", new eckit::OStreamTarget(eckit::Log::info()))));
  } else {
    testChannel_.reset(new eckit::Channel());
  }
  return *testChannel_;
}

eckit::Channel& LibOOPS::infoChannel() const {
  if (rank_ == 0) {return eckit::Log::info();}
  if (!infoChannel_) infoChannel_.reset(new eckit::Channel());
  return *infoChannel_;
}

eckit::Channel& LibOOPS::debugChannel() const {
  if (debugChannel_) {return *debugChannel_;}
  if (debug_) {
    debugChannel_.reset(new eckit::Channel(new eckit::PrefixTarget(predebug_)));
  } else {
    debugChannel_.reset(new eckit::Channel());
  }
  return *debugChannel_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

