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

#include <iostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/LibOOPS.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCountHelper.h"
#include "oops/util/printRunStats.h"
#include "oops/util/Stacktrace.h"
#include "oops/util/TimerHelper.h"
#include "oops/util/workflow.h"

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

  // Command line options and arguments
  std::string infilename = "";
  std::string outfilename = "";
  for (int i = 1; i < argc; ++i) {
    std::string item = static_cast<std::string>(argv[i]);
    if (item == "-h" || item == "--help") {
      eckit::PathName argv0 = static_cast<std::string>(argv[0]);
      argv0 = argv0.baseName();
      std::cout << "Usages:" << std::endl;
      std::cout << "  # run main application:" << std::endl;
      std::cout << "  " << argv0 << " input-file [output-file]" << std::endl;
      std::cout << "  # run main application without validating YAML file:" << std::endl;
      std::cout << "  " << argv0 << " --no-validate input-file" << std::endl;
      std::cout << "  # check input YAML file against its schema:" << std::endl;
      std::cout << "  " << argv0 << " --validate-only input-file" << std::endl;
      std::cout << "  # write input file schema to given file name:" << std::endl;
      std::cout << "  " << argv0 << " --output-json-schema=file-name" << std::endl;
      std::cout << "  # print this help and exit:" << std::endl;
      std::cout << "  " << argv0 << " --help" << std::endl;
      is_print_help_only_ = true;
    } else if (item == "--validate-only") {
      is_validate_only_ = true;
    } else if (item == "--no-validate") {
      validate_ = false;
    } else if (item.rfind("--output-json-schema=", 0) == 0) {
      output_json_schema_path_ = item.substr(item.find("=") + 1);
    } else if (infilename.empty()) {
      infilename = item;
    } else if (outfilename.empty()) {
      outfilename = item;
    } else {
      ABORT(item + ": unknown option or positional argument");
    }
  }

  if (infilename.empty()) {
    if (!is_print_help_only_ && output_json_schema_path_.empty()) {
      ABORT("Positional argument 1 must be the file name of a YAML configuration file");
    }
  } else if (!output_json_schema_path_.empty()) {
    ABORT("Cannot specify config file name (positional argument 1) with --output-json-schema=...");
  } else if (is_validate_only_ && !outfilename.empty()) {
    ABORT("Cannot specify output file name (positional argument 2) with --validate-only");
  } else {
    // Read configuration
    eckit::PathName infilepathname = infilename;
    config_.reset(new eckit::YAMLConfiguration(infilepathname));

    if (!is_validate_only_) {
      // Get configuration file and optional output file from command line
      LibOOPS::instance().initialise();

      // Configure TestReference with "test:" sub-config
      if (config_->has("test"))
        LibOOPS::instance().testReferenceInitialise(config_->getSubConfiguration("test"));

      if (!outfilename.empty()) {
        eckit::PathName outfilepathname = outfilename;
        LibOOPS::instance().teeOutput(outfilepathname);
      }
    }

    Log::info() << "Configuration input file is: " << infilepathname << std::endl;
    Log::info() << "Full configuration is:"  << *config_ << std::endl;
  }
}

// -----------------------------------------------------------------------------

Run::~Run() {
  if (!is_print_help_only_ && !is_validate_only_ && output_json_schema_path_.empty()) {
    LibOOPS::instance().finalise();  // Finalize MPI and logs
  }
}

// -----------------------------------------------------------------------------

int Run::execute(const Application & app, const eckit::mpi::Comm & comm) {
  if (is_print_help_only_) {
    return 0;
  }
  int status = 1;
  try {
    if (!output_json_schema_path_.empty()) {
      app.outputSchema(output_json_schema_path_);
      Log::info() << "Output JSON Schema file: " << output_json_schema_path_ << std::endl;
      status = 0;
    } else if (is_validate_only_) {
      app.validateConfig(*config_);
      Log::info() << "Configuration OK" << std::endl;
      status = 0;
    } else {
      // Start measuring performance
      util::TimerHelper::setComm(comm);
      util::TimerHelper::start();
      util::ObjectCountHelper::start();
      util::printRunStats("Run start", true, comm);
      if (config_->getBool("ecflow", false)) util::use_ecflow();
      // Run application
      Log::info() << "Run: Starting " << app << std::endl;
      status = app.execute(*config_, validate_);
      Log::info() << std::endl << "Run: Finishing " << app << std::endl;
      // Performance diagnostics
      util::ObjectCountHelper::stop();
      util::TimerHelper::stop();
      util::printRunStats("Run end", true, comm);
      Log::info() << "Run: Finishing " << app << " with status = " << status << std::endl;
    }
  }
  catch(const eckit::Exception & e) {
    status = 1;
    Log::error() << e.what() << " caught in "  << Here() << std::endl;
    Log::error() << "Exception: " << app << " terminating..." << std::endl;
    eckit::Exception::exceptionStack(eckit::Log::error(), true);
    util::unwind_exception_stack(e, eckit::Log::error());
  }
  catch(const std::exception & e) {
    status = 1;
    Log::error() << "Exception: " << e.what() << std::endl;
    Log::error() << "Exception: " << app << " terminating..." << std::endl;
    util::unwind_exception_stack(e, eckit::Log::error());
  }
  catch(...) {
    status = 1;
    Log::error() << "Unknown exception: " << app << " terminating..." << std::endl;
  }
  return status;
}

// -----------------------------------------------------------------------------

}  // namespace oops
