/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_RUN_H_
#define OOPS_RUNS_RUN_H_

#include <boost/scoped_ptr.hpp>

#include "eckit/runtime/Main.h"
#include "eckit/config/YAMLConfiguration.h"
#include "util/Timer.h"

namespace oops {
  class Application;

/*!
 *  Run encapsulates one OOPS run.
 */

// -----------------------------------------------------------------------------

class Run : public eckit::Main {
 public:
  Run(int argc, char** argv);
  virtual ~Run();
  void execute(const Application &);

 protected:
  const eckit::Configuration & config() const {return *config_;}

 private:
  boost::scoped_ptr<const eckit::YAMLConfiguration> config_;
  boost::scoped_ptr<util::Timer> timer_;
};

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_RUNS_RUN_H_
