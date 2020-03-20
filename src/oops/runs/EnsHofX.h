/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSHOFX_H_
#define OOPS_RUNS_ENSHOFX_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/HofX.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class EnsHofX : public Application {
 public:
// -----------------------------------------------------------------------------
  explicit EnsHofX(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateObsFilterFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~EnsHofX() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup initial states
    const eckit::LocalConfiguration initialConfig(fullConfig, "Initial Condition");
    std::vector<eckit::LocalConfiguration> memberConf;
    initialConfig.get("state", memberConf);
    Log::debug() << "EnsHofX: using " << memberConf.size() << " states." << std::endl;

//  Get the MPI partition information
    const int ntasks = oops::mpi::comm().size();
    const int mytask = oops::mpi::comm().rank();
    int members = 1;
    int tasks_per_member = 0;
    int mymember = 0;

    if ( fullConfig.has("EnsembleApplication.members") &&
         !(fullConfig.has("EnsembleApplication.current_member")) ) {
      members = fullConfig.getInt("EnsembleApplication.members");
      tasks_per_member = ntasks / members;
      mymember = mytask / tasks_per_member + 1;
      Log::info() << "Running " << members << " EnsForecast members handled by "
                  << ntasks << " total MPI tasks and "
                  << tasks_per_member << " MPI tasks per member." << std::endl;
    } else if ( !(fullConfig.has("EnsembleApplication.members")) &&
               fullConfig.has("EnsembleApplication.current_member") ) {
      tasks_per_member = ntasks;
      mymember = fullConfig.getInt("EnsembleApplication.current_member");
      Log::info() << "Running EnsForecast member number " << mymember
                  << " handled by " << ntasks << " total MPI tasks." << std::endl;
    } else {
      ABORT("The options are EnsembleApplication.current_member OR EnsembleApplication.members");
    }

    ASSERT(ntasks%members == 0);

//  Create  the communicator for each member, named comm_member_{i}
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = eckit::mpi::comm().split(mymember, commName);

  // Add the useful info in the eckit configuration
    eckit::LocalConfiguration config(fullConfig);
    config.set("Observations.member", mymember);
    config.set("Initial Condition", memberConf[mymember-1]);

    Log::debug() << "EnsHofX config for member 0 = " << config << std::endl;

    HofX<MODEL> hofx(commMember);
    hofx.execute(config);

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EnsHofX<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_ENSHOFX_H_
