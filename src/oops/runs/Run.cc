/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/runs/Run.h"

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/runs/Application.h"
#include "oops/util/LibOOPS.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCountHelper.h"
#include "oops/util/TimerHelper.h"

#ifdef ENABLE_GPTL
#include <gptl.h>
#endif

// GPTL implementation in JEDI requires retrieval of both integer and string env. var.
int getEnv(const std::string& env, int default_value) {
  if (::getenv(env.c_str())) {return eckit::Translator<std::string, int>()(::getenv(env.c_str()));}
  return default_value;
}

std::string getEnv(const std::string& env, std::string default_value) {
  if ((::getenv(env.c_str()))) {return ::getenv(env.c_str());}
  return default_value;
}

namespace oops {

  // -----------------------------------------------------------------------------

Run::Run(int argc, char** argv) : eckit::Main(argc, argv, "OOPS_HOME"), config_(), timer_() {
  // Initialize MPI and LibOOPS variables that require eckit::Main
#ifdef ENABLE_GPTL
  int do_profile = getEnv("OOPS_PROFILE", 0);  // Default is profiling disabled
  if (do_profile) {
    int ret;
    // All GPTLsetutr() and GPTLsetoption() calls must precede GPTLinitialize()

    // Default JEDI behavior: Change the GPTL default from expensive gettimeofday() to inexpensive
    // register read ("nanotime"). User can control with $OOPS_UNDERLYING_TIMER
    std::string utr = getEnv("OOPS_UNDERLYING_TIMER", "nanotime");
    if (utr == "nanotime") {
      ret = GPTLsetutr(GPTLnanotime);      // Use x86-specific register read to gather timings
    } else if (utr == "gettimeofday") {
      ret = GPTLsetutr(GPTLgettimeofday);  // Use gettimeofday() to gather timings
    } else {
      Log::warning() << "OOPS_UNDERLYING_TIMER=" << utr << " is invalid: ignoring" << std::endl;
    }

    // Set printout method. Another common method is full_tree, but that can create tons of output
    ret = GPTLsetoption(GPTLprint_method, GPTLmost_frequent);

    // For MPI codes, synchronize and time collectives and receives to avoid mis-assignment
    // of load imbalance to true MPI time. Disable by changing 1 to 0 here or via env. var.
    int sync_mpi = getEnv("OOPS_SYNC_MPI", 1);
    if (sync_mpi)
      ret = GPTLsetoption(GPTLsync_mpi, 1);

    // Default behavior is do not print hi-water mark increases as the process runs
    int dopr_memusage = getEnv("OOPS_MEMUSAGE", 0);
    if (dopr_memusage)
      ret = GPTLsetoption(GPTLdopr_memusage, 1);   // Print growth of resident set size (RSS)

    ret = GPTLinitialize();  // Initialize GPTL timing library
  }
#endif

  LibOOPS::instance().initialise();

// Get configuration file and optional output file from command line
  ASSERT(argc >= 2);
  eckit::PathName configfile = argv[1];
  if (argc == 3) {
    eckit::PathName outputfile;
    outputfile = argv[2];
    LibOOPS::instance().teeOutput(outputfile);
  }

// Read configuration
  config_.reset(new eckit::YAMLConfiguration(configfile));


// Read tolerances for testing
  if (config_->has("test")) {
    auto testConf = config_->getSubConfiguration("test");
    LibOOPS::instance().testReferenceInitialise(testConf);
  }

  Log::info() << "Configuration input file is: " << configfile << std::endl;
  Log::info() << "Full configuration is:"  << *config_ << std::endl;

// Start measuring performance
  util::ObjectCountHelper::start();
  util::TimerHelper::start();
}

// -----------------------------------------------------------------------------

Run::~Run() {
  LibOOPS::instance().finalise();  // Finalize MPI and logs
}

// -----------------------------------------------------------------------------

int Run::execute(const Application & app) {
  int status = 1;
  Log::info() << "Run: Starting " << app << std::endl;
  try {
    status = app.execute(*config_);
  }
  catch(const eckit::Exception & e) {
    status = 1;
    Log::error() << e.what() << " caught in "  << Here() << std::endl;
    Log::error() << "Exception: " << app << " terminating..." << std::endl;
    eckit::Exception::exceptionStack(eckit::Log::error(), true);
  }
  catch(const std::exception & e) {
    status = 1;
    Log::error() << "Exception: " << e.what() << std::endl;
    Log::error() << "Exception: " << app << " terminating..." << std::endl;
  }
  catch(...) {
    status = 1;
    Log::error() << "Unknown exception: " << app << " terminating..." << std::endl;
  }
  Log::info() << "Run: Finishing " << app << std::endl;

// Performance diagnostics
  util::ObjectCountHelper::stop();
  util::TimerHelper::stop();

  Log::info() << "Run: Finishing " << app << " with status = " << status << std::endl;
  return status;
}

// -----------------------------------------------------------------------------

}  // namespace oops
