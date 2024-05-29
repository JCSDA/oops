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
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/log/Log.h"
#include "eckit/log/OStreamTarget.h"
#include "eckit/log/PrefixTarget.h"
#include "eckit/utils/Translator.h"

#include "oops/mpi/mpi.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/LibOOPS.h"
#include "oops/util/Logger.h"
#include "oops/util/signal_trap.h"

#ifdef ENABLE_GPTL
#include <gptl.h>
#include <gptlmpi.h>

int do_profile = false;  // Flag says whether to enable profiling with GPTL (default false)
#endif

namespace oops {

//------------------------------------------------------------------------------

  bool getEnv(const std::string& env, bool default_value) {
    if (::getenv(env.c_str()))
      {return eckit::Translator<std::string, bool>()(::getenv(env.c_str()));}
    return default_value;
  }

  int getEnv(const std::string& env, int default_value) {
    if (::getenv(env.c_str()))
      {return eckit::Translator<std::string, int>()(::getenv(env.c_str()));}
    return default_value;
  }

  std::vector<int> getEnvList(const std::string& env, std::vector<int> default_value) {
    if (::getenv(env.c_str())) {
      const std::string pelist(::getenv(env.c_str()));
      const std::set<int> pes = parseIntSet(pelist);
      std::vector<int> info_pes;
      std::copy(pes.begin(), pes.end(), std::back_inserter(info_pes));
      return info_pes;
    }
    return default_value;
  }

//------------------------------------------------------------------------------

static LibOOPS liboops;

LibOOPS::LibOOPS() : Library("oops"),
                     initialised_(false), rank_(0), test_(false),
                     info_(false), preinfo_(""),
                     debug_(false), predebug_("OOPS_DEBUG"),
                     trace_(false), pretrace_("OOPS_TRACE"),
                     prestat_("OOPS_STATS") {
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
  std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char nowstr[100];
  std::strftime(nowstr, sizeof(nowstr), "%F %T (UTC%z)", std::localtime(&now));
  rank_ = oops::mpi::world().rank();
  if (rank_ == 0) Log::info() << "OOPS Starting " << nowstr << std::endl;

  std::vector<int> iis(getEnvList("OOPS_INFO", std::vector<int>(1, 0)));
  if (rank_ == 0) info_ = true;  // INFO is always logged for rank 0
  for (size_t iter = 0; iter < iis.size(); iter++) {
    if (iis[iter] < 0) info_ = true;
    if (static_cast<int>(rank_) == iis[iter]) info_ = true;
  }
  if (info_ && rank_ > 0) {
    preinfo_ += "[" + std::to_string(rank_) + "]";
    prestat_ += "[" + std::to_string(rank_) + "]";
  }

  const int it = getEnv("OOPS_TRACE", 0);
  if (it > 0) trace_ = info_;
  if (it < 0) trace_ = true;
  pretrace_ += "[" + std::to_string(rank_) + "]";

  const int id = getEnv("OOPS_DEBUG", 0);
  if (id > 0) debug_ = info_;
  if (id < 0) debug_ = true;
  predebug_ += "[" + std::to_string(rank_) + "]";

  const int do_trapfpe = getEnv("OOPS_TRAPFPE", 0);
  int do_abortfpe;
  if (do_trapfpe) {
    // If SIGFPE trapping is enabled, default is to abort. Use caution trapping but not aborting:
    // It can possibly result in gargantuan output to stderr
    do_abortfpe = getEnv("OOPS_ABORTFPE", 1);
    util::trap_sigfpe(do_abortfpe);
  }

#ifdef ENABLE_GPTL
  do_profile = getEnv("OOPS_PROFILE", 0);
#endif

  initialised_ = true;
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

void LibOOPS::testReferenceInitialise(const eckit::LocalConfiguration &testConf) {
  // testStream_ is used by TestReference for comparing test output
  // with a reference file
  if ( rank_ == 0 ) {
    test_ = true;
    testChannel().addStream(testStream_);
    testReference_.initialise(testConf);
  }
}

/** Clears logs and finalises MPI (unless \p finaliseMPI is false).
 * To be called in on leaving `main()` by the destructor of `oops::Run`.
 */
void LibOOPS::finalise(bool finaliseMPI) {
    eckit::Log::flush();
#ifdef ENABLE_GPTL
    if (do_profile) {
      int ret;
      if (oops::mpi::world().rank() == 0) {
        oops::Log::info() << "Calling GPTLpr(0)" << std::endl;
        ret = GPTLpr(0);                       // Print timing info for rank 0
      }
      // Summarize timing info across ranksL For now assume MPI_COMM_WORLD
      if (oops::mpi::world().size() > 1) {
        oops::Log::info() << "Calling GPTLpr_summary" << std::endl;
        ret = GPTLpr_summary(MPI_COMM_WORLD);
      }
    }
#endif

    if ( rank_ == 0 ) {
      testReference_.finalise(testStream_.str());
    }

    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char nowstr[100];
    std::strftime(nowstr, sizeof(nowstr), "%F %T (UTC%z)", std::localtime(&now));
    Log::info() << "OOPS Ending   " << nowstr << std::endl;

    // Make sure that these specialised channels that wrap eckit::Log::info() are
    // destroyed before eckit::Log::info gets destroyed.
    // Just in case someone still tries to log, we reset to empty channels.

    debugChannel_.reset(new eckit::Channel());
    traceChannel_.reset(new eckit::Channel());
    statsChannel_.reset(new eckit::Channel());
    testChannel_.reset(new eckit::Channel());
    // Destroy info channel last after other channels have flushed all output
    infoChannel_.reset(new eckit::Channel());

    if (finaliseMPI) eckit::mpi::finaliseAllComms();
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
  // use pre-main channel if library was not initialized
  if (!initialised_) {return eckit::Log::info();}
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
  if (rank_ == 0 || info_) {
    statsChannel_.reset(new eckit::Channel(
      new eckit::PrefixTarget(prestat_, new eckit::OStreamTarget(eckit::Log::info()))));
  } else {
    statsChannel_.reset(new eckit::Channel());
  }
  return *statsChannel_;
}

eckit::Channel& LibOOPS::testChannel() const {
  // use pre-main channel if library was not initialized
  if (!initialised_) {return eckit::Log::info();}
  if (testChannel_) {return *testChannel_;}
  if (test_) {
    testChannel_.reset(new eckit::Channel(
      new eckit::PrefixTarget("Test     :", new eckit::OStreamTarget(eckit::Log::info()))));
  } else {
    testChannel_.reset(new eckit::Channel());
  }
  testChannel_->setf(std::ios::scientific);
  testChannel_->precision(std::numeric_limits<double>::digits10+1);
  return *testChannel_;
}

eckit::Channel& LibOOPS::infoChannel() const {
  // use pre-main channel if library was not initialized
  if (!initialised_) {return eckit::Log::info();}
  if (infoChannel_) {return *infoChannel_;}
  if (rank_ == 0) {
    return eckit::Log::info();
  } else if (info_ && rank_ > 0) {
    infoChannel_.reset(new eckit::Channel(
      new eckit::PrefixTarget(preinfo_, new eckit::OStreamTarget(eckit::Log::info()))));
  } else {
    infoChannel_.reset(new eckit::Channel());
  }
  if (!infoChannel_) infoChannel_.reset(new eckit::Channel());
  return *infoChannel_;
}

eckit::Channel& LibOOPS::debugChannel() const {
  // use pre-main channel if library was not initialized
  if (!initialised_) {return eckit::Log::info();}
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
