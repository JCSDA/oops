/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_EDA_H_
#define OOPS_RUNS_EDA_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/Variational.h"
#include "oops/util/Logger.h"

namespace oops {

template<typename MODEL, typename OBS> class EDA : public Application {
 public:
// -----------------------------------------------------------------------------
  explicit EDA(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~EDA() {}
// -----------------------------------------------------------------------------
int execute(const eckit::Configuration & fullConfig) const {
// Get the MPI partition information
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
    Log::info() << "Running " << members << " EDA members handled by "
                << ntasks << " total MPI tasks and "
                << tasks_per_member << " MPI tasks per member." << std::endl;
  } else if ( !(fullConfig.has("EnsembleApplication.members")) &&
             fullConfig.has("EnsembleApplication.current_member") ) {
    tasks_per_member = ntasks;
    mymember = fullConfig.getInt("EnsembleApplication.current_member");
    Log::info() << "Running EDA member number " << mymember
                << " handled by " << ntasks << " total MPI tasks." << std::endl;
  } else {
    ABORT("The options are EnsembleApplication.current_member OR EnsembleApplication.members");
  }

  ASSERT(ntasks%members == 0);

// Create  the communicator for each member, named comm_member_{i}
  std::string commNameStr = "comm_member_" + std::to_string(mymember);
  char const *commName = commNameStr.c_str();
  eckit::mpi::Comm & commMember = eckit::mpi::comm().split(mymember, commName);

// Add the useful info in the eckit configuration
  eckit::LocalConfiguration config(fullConfig);

// To write analysis in different files
  config.set("output.member", mymember);

// To read different backgrounds
  config.set("cost_function.Jb.Background.member", mymember);

// To perturb the observations and write hofx in different files:
  config.set("cost_function.Jo.member", mymember);
  if (mymember > 1) config.set("cost_function.Jo.ObsPert", true);

  Log::debug() << "EDA config for member 0 = " << config << std::endl;

  Variational<MODEL, OBS> var(commMember);
  return var.execute(config);
}
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EDA<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_EDA_H_
