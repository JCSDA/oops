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

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the EnsembleApplication application.
template <typename APP>
class EnsembleApplicationParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(EnsembleApplicationParameters, ApplicationParameters)

 public:
  /// Parameters containing a list of YAML files for each ensemble member to be processed.
  RequiredParameter<std::vector<std::string>> files{"files", this};
};

// -----------------------------------------------------------------------------

template <typename APP>
class EnsembleApplication : public Application {
  typedef EnsembleApplicationParameters<APP> EnsembleApplicationParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit EnsembleApplication(const eckit::mpi::Comm & comm = oops::mpi::world()) :
      Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~EnsembleApplication() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    EnsembleApplicationParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Get the list of YAML files
    const std::vector<std::string> &files = params.files.value();

    Log::info() << "EnsembleApplication YAML files:" << files << std::endl;

//  Get the MPI partition
    const int nmembers = files.size();
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

//  Each member uses a different configuration:
    eckit::PathName confPath = files[mymember-1];
    eckit::YAMLConfiguration memberConf(confPath);

    Log::debug() << "EnsembleApplication config for member " << mymember << ": "
                 << memberConf << std::endl;

    APP ensapp(commMember);
    return ensapp.execute(memberConf);
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    EnsembleApplicationParameters_ params;
    params.outputSchema(outputPath);
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
