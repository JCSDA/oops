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

#include <boost/scoped_ptr.hpp>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/runs/Application.h"
#include "util/ObjectCountHelper.h"
#include "util/TimerHelper.h"
#include "util/Logger.h"
#include "util/LibOOPS.h"

extern "C" {
    void mpl_start_f90();
    void mpl_end_f90();
}

namespace oops {

// -----------------------------------------------------------------------------

Run::Run(int argc, char** argv) : eckit::Main(argc, argv, "OOPS_HOME"), config_(), timer_() {
// Run MPI for NICAS
  mpl_start_f90();

// Get configuration file from command line
  ASSERT(argc >= 2);
  eckit::PathName configfile = argv[argc - 1];

// Read configuration
  config_.reset(new eckit::YAMLConfiguration(configfile));

  Log::info() << "Configuration input file is: " << configfile << std::endl;
  Log::info() << "Full configuration is:"  << *config_ << std::endl;

// Start measuring performance
  util::ObjectCountHelper::start();
  util::TimerHelper::start();
}

// -----------------------------------------------------------------------------

Run::~Run() {
// Finalize MPI for NICAS
    mpl_end_f90();

    LibOOPS::instance().finalise();
}

// -----------------------------------------------------------------------------

void Run::execute(const Application & app) {
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
  if (status) ::exit(status);
}

// -----------------------------------------------------------------------------

}  // namespace oops
