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
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/Forecast.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class EnsForecast : public Application {
  typedef Geometry<MODEL>         Geometry_;
  typedef Model<MODEL>            Model_;
  typedef ModelAuxControl<MODEL>  ModelAux_;
  typedef State<MODEL>            State_;

 public:
// -----------------------------------------------------------------------------
  explicit EnsForecast(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~EnsForecast() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    std::vector<eckit::LocalConfiguration> memberConf;
    fullConfig.get("members", memberConf);

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

    Log::info() << members << " EnsForecast members handled by " << ntasks
                << " total MPI tasks and " << tasks_per_member << " MPI tasks per member"
                << std::endl;

// Add the useful info in the eckit configuration
    eckit::LocalConfiguration config(fullConfig);
    config.set("output.member", mymember);  // To write analysis in different files
    config.set("initial", memberConf[mymember-1]);

    Log::debug() << "EnsForecast config for member 0 = " << config << std::endl;

    Forecast<MODEL> fc(commMember);
    fc.execute(config);

    return 0;
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
