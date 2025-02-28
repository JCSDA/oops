/*
 * (C) Copyright 2018-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSEMBLEAPPLICATION_H_
#define OOPS_RUNS_ENSEMBLEAPPLICATION_H_

#include <string>
#include <vector>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename APP>
class EnsembleApplication : public Application {
 public:
// -----------------------------------------------------------------------------
  explicit EnsembleApplication(const eckit::mpi::Comm & comm = oops::mpi::world()) :
      Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~EnsembleApplication() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Assert that the config does not contain both files and members
    if (fullConfig.has("files") && fullConfig.has("members")) {
      throw eckit::BadParameter("EnsembleApplication: 'files' and 'members' "
                                "cannot be specified at the same time", Here());
    }

//  Get the list of configurations
    std::vector<eckit::LocalConfiguration> memberConfigs;
    if (fullConfig.has("files")) {
      // Read the configurations from the files
      for (const std::string &file : fullConfig.getStringVector("files")) {
        const eckit::PathName yamlPathFile(file);
        const eckit::YAMLConfiguration memberConfig(yamlPathFile);
        memberConfigs.push_back(eckit::LocalConfiguration(memberConfig));
      }
    } else if (fullConfig.has("members")) {
      // Copy the configurations from the members
      memberConfigs = fullConfig.getSubConfigurations("members");
    } else {
      throw eckit::BadParameter("EnsembleApplication: either 'files' or 'members' "
                                "must be specified in the configuration", Here());
    }

//  Get the MPI partition
    const int nmembers = memberConfigs.size();
    const int ntasks = this->getComm().size();
    const int mytask = this->getComm().rank();
    const int tasks_per_member = ntasks / nmembers;
    const int mymember = mytask / tasks_per_member + 1;

    Log::info() << "Running " << nmembers << " EnsembleApplication members handled by "
                << ntasks << " total MPI tasks and "
                << tasks_per_member << " MPI tasks per member." << std::endl;

    ASSERT(ntasks%nmembers == 0);

//  Create the communicator for each member, named comm_member_{i}:
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = this->getComm().split(mymember, commName);

//  Run each member with the corresponding configuration:
    APP ensapp(commMember);
    return ensapp.execute(memberConfigs[mymember-1]);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsembleApplication<>";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ENSEMBLEAPPLICATION_H_
