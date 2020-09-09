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
#include <iomanip>
#include <sstream>
#include <string>

#include "eckit/log/Log.h"
#include "eckit/log/OStreamTarget.h"
#include "eckit/log/PrefixTarget.h"
#include "eckit/utils/Translator.h"

#include "oops/mpi/mpi.h"
#include "oops/util/LibOOPS.h"

extern void trap_sigfpe(const int);

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

LibOOPS::LibOOPS() : Library("oops"), rank_(0),
                     debug_(false), predebug_("OOPS_DEBUG"),
                     trace_(false), pretrace_("OOPS_TRACE") {
}

LibOOPS::~LibOOPS() {
}

LibOOPS& LibOOPS::instance() {
  return liboops;
}

/** Initialization of MPI and dependent variables.
 * To be called in `main()` by constructor of `oops::Run`.  This method initializes MPI and 
 * associated variables that must be initialized after static-init time, and only once `eckit::Main`
 * has been created.
 */
void LibOOPS::initialise() {
  rank_ = oops::mpi::world().rank();

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
  const int do_trapfpe = getEnv("OOPS_TRAPFPE", 0);
  int do_abortfpe;
  if (do_trapfpe) {
    // If SIGFPE trapping is enabled, default is to abort. Use caution trapping but not aborting:
    // It can possibly result in gargantuan output to stderr
    do_abortfpe = getEnv("OOPS_ABORTFPE", 1);
    trap_sigfpe(do_abortfpe);
  }
}

/** Add a rank-dependent tee file
 * To be called in `main()` by constructor of `oops::Run` when the outputfile argument is specified
 */
void LibOOPS::teeOutput(const std::string & fileprefix) {
  std::string teefile = fileprefix;
  if (rank_ != 0) {
    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << rank_;
    teefile = teefile + "." + ss.str();
  }
  eckit::Log::addFile(teefile);
}

/** Clears logs and finalises MPI (unless \p finaliseMPI is false).
 * To be called in on leaving `main()` by the destructor of `oops::Run`.
 */
void LibOOPS::finalise(bool finaliseMPI) {
    eckit::Log::flush();

    // Make sure that these specialised channels that wrap eckit::Log::info() are
    // destroyed before eckit::Log::info gets destroyed.
    // Just in case someone still tries to log, we reset to empty channels.
    infoChannel_.reset(new eckit::Channel());
    debugChannel_.reset(new eckit::Channel());
    traceChannel_.reset(new eckit::Channel());
    statsChannel_.reset(new eckit::Channel());
    testChannel_. reset(new eckit::Channel());

    if (finaliseMPI)
      eckit::mpi::finaliseAllComms();
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

