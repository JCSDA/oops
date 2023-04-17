/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/coupled/UtilsCoupled.h"
#include "eckit/exception/Exceptions.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
std::vector<Variables> splitVariables(const Variables & varsToSplit,
                                      const std::vector<Variables> & varsOfCoupledModels) {
  std::vector<Variables> splitvars = varsOfCoupledModels;
  if (varsOfCoupledModels.size() > 0) {
    // check that the same variable isn't specified in both models' variables
    Variables allvars;
    size_t allvarsize = 0;
    for (const auto & vars : varsOfCoupledModels) {
      allvars += vars;
      allvarsize += vars.size();
    }
    if (allvars.size() != allvarsize) {
      Log::error() << "Variables for different components of coupled model can not overlap."
                   << std::endl;
      throw eckit::BadParameter("Variables for different components of coupled "
                                "model can not overlap", Here());
    }
    // decide what variables are provided by what model
    for (auto & element : splitvars) {
      element.intersection(varsToSplit);
    }
    // check that all variables are accounted for
    Variables allreqvars;
    for (const auto & element : splitvars) {
      allreqvars += element;
    }
    if (allreqvars != varsToSplit) {
      Log::error() << "Not all variables can be provided by the coupled model. " << std::endl
                   << "Requested: " << varsToSplit << std::endl
                   << "Available: " << allreqvars << std::endl;
      throw eckit::UserError("Not all variables can be provided by the coupled model");
    }
  }
  return splitvars;
}

}  // namespace oops
