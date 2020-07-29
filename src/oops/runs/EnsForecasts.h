/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_ENSFORECASTS_H_
#define OOPS_RUNS_ENSFORECASTS_H_

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/Forecast.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class EnsForecast : public Application {
 public:
// -----------------------------------------------------------------------------
  explicit EnsForecast(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~EnsForecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    std::vector<eckit::LocalConfiguration> memberConf;
    fullConfig.get("members", memberConf);

// Get the MPI partition information
    const int ntasks = oops::mpi::comm().size();
    const int mytask = oops::mpi::comm().rank();
    int members = 1;
    int tasks_per_member = 0;
    int mymember = 0;

    if ( fullConfig.has("ensemble application.members") &&
         !(fullConfig.has("ensemble application.current member")) ) {
      members = fullConfig.getInt("ensemble application.members");
      tasks_per_member = ntasks / members;
      mymember = mytask / tasks_per_member + 1;
      Log::info() << "Running " << members << " EnsForecast members handled by "
                  << ntasks << " total MPI tasks and "
                  << tasks_per_member << " MPI tasks per member." << std::endl;
    } else if ( !(fullConfig.has("ensemble application.members")) &&
               fullConfig.has("ensemble application.current member") ) {
      tasks_per_member = ntasks;
      mymember = fullConfig.getInt("ensemble application.current member");
      Log::info() << "Running EnsForecast member number " << mymember
                  << " handled by " << ntasks << " total MPI tasks." << std::endl;
    } else {
      ABORT("The options are ensemble application.current member OR ensemble application.members");
    }

    ASSERT(ntasks%members == 0);

// Create  the communicator for each member, named comm_member_{i}
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = eckit::mpi::comm().split(mymember, commName);

// Add the useful info in the eckit configuration
    eckit::LocalConfiguration config(fullConfig);
    config.set("output.member", mymember);  // To write analysis in different files
    config.set("initial condition", memberConf[mymember-1]);

    Log::debug() << "EnsForecast config for member 0 = " << config << std::endl;

    Forecast<MODEL> fc(commMember);
    return fc.execute(config);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EnsForecast<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ENSFORECASTS_H_
