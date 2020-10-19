/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSEMBLEAPPLICATION_H_
#define OOPS_RUNS_ENSEMBLEAPPLICATION_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename APP>
class EnsembleApplication : public Application {
 public:
// -----------------------------------------------------------------------------
  explicit EnsembleApplication(const eckit::mpi::Comm & comm = oops::mpi::world()) :
      Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~EnsembleApplication() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
  // Get the list of yaml files
    std::vector<std::string> listConf;
    fullConfig.get("files", listConf);
    Log::info() << "EnsembleApplication yaml files:" << listConf << std::endl;

  // Get the MPI partition
    const int nmembers = listConf.size();
    const int ntasks = this->getComm().size();
    const int mytask = this->getComm().rank();
    const int tasks_per_member = ntasks / nmembers;
    const int mymember = mytask / tasks_per_member + 1;

    Log::info() << "Running " << nmembers << " EnsembleApplication members handled by "
                << ntasks << " total MPI tasks and "
                << tasks_per_member << " MPI tasks per member." << std::endl;

    ASSERT(ntasks%nmembers == 0);

  // Create  the communicator for each member, named comm_member_{i}:
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = this->getComm().split(mymember, commName);

  // Each member uses a different configuration:
    eckit::PathName confPath = listConf[mymember-1];
    eckit::YAMLConfiguration memberConf(confPath);

    Log::debug() << "EnsembleApplication config for member " << mymember << ": "
                 << memberConf << std::endl;

    APP ensapp(commMember);
    return ensapp.execute(memberConf);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EnsembleApplication<>";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ENSEMBLEAPPLICATION_H_
