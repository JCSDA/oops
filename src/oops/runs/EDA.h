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

template<typename MODEL> class EDA : public Application {
 public:
// -----------------------------------------------------------------------------
  explicit EDA(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~EDA() {}
// -----------------------------------------------------------------------------
int execute(const eckit::Configuration & fullConfig) const {
// Create  the communicator for each member, named comm_member_{i}
  const int members = fullConfig.getInt("EnsembleApplication.members");
  const int ntasks = oops::mpi::comm().size();
  const int mytask = oops::mpi::comm().rank();
  const int tasks_per_member = ntasks / members;
  const int mymember = mytask / tasks_per_member + 1;
  ASSERT(ntasks%members == 0);

  std::string commNameStr = "comm_member_" + std::to_string(mymember);
  char const *commName = commNameStr.c_str();
  eckit::mpi::Comm & commMember = eckit::mpi::comm().split(mymember, commName);

  Log::info() << members << " EDA members handled by " << ntasks << " total MPI tasks and "
              << tasks_per_member << " MPI tasks per member" << std::endl;

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

  Variational<MODEL> var(commMember);
  var.execute(config);

  return 0;
}
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EDA<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_EDA_H_
