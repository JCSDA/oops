/*
 * (C) Copyright 2021- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_EMPTYPARAMETERS_H_
#define OOPS_UTIL_PARAMETERS_EMPTYPARAMETERS_H_

#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief A subclass of Parameters storing no options.
/// It can be used when a Parameters class needs to be defined, and no options
/// are required.
class EmptyParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(EmptyParameters, Parameters)
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_EMPTYPARAMETERS_H_
