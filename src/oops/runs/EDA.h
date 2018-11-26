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
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/Variational.h"
#include "oops/util/Logger.h"

namespace oops {

template<typename MODEL> class EDA : public Application {
 public:
// -----------------------------------------------------------------------------
  EDA() {}
// -----------------------------------------------------------------------------
  virtual ~EDA() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    const int members = fullConfig.getInt("EDA.members");
    const int ntasks = oops::mpi::comm().size();
    const int mytask = oops::mpi::comm().rank();
    Log::debug() << "EDA members = " << members << ", ntasks = " << ntasks
                 << ", mytask = " << mytask << std::endl;

//  Should create members groups but for now just one task per member
    ASSERT(ntasks == members);

    eckit::LocalConfiguration config(fullConfig);
    config.set("cost_function.Jo.member", mytask);
    if (mytask > 0) config.set("cost_function.Jo.ObsPert", true);
    Log::debug() << "EDA config = " << config << std::endl;

    Variational<MODEL> var;
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
