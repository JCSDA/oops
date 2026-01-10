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

#include <memory>
#include <string>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/runtime/Main.h"
#include "oops/mpi/mpi.h"

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
  int execute(const Application &, const eckit::mpi::Comm & comm = oops::mpi::world());

 protected:
  const eckit::Configuration & config() const {return *config_;}

 private:
  bool is_print_help_only_ = false;
  std::unique_ptr<const eckit::YAMLConfiguration> config_;
};

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_RUNS_RUN_H_
