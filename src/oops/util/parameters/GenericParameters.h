/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_GENERICPARAMETERS_H_
#define OOPS_UTIL_PARAMETERS_GENERICPARAMETERS_H_

#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief A subclass of Parameters storing the values of all options in an
/// eckit::LocalConfiguration object.
///
/// This object can be accessed by calling the `value()` method of the `config` member variable.
///
/// Values loaded into a ConfigurationParameter aren't validated in any way; parts of the code
/// using GenericParameters should therefore ideally be refactored, replacing this class with a
/// dedicated subclass of Parameters storing each parameter in a separate
/// (Optional/Required)Parameter object.
class GenericParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(GenericParameters, Parameters)

 public:
  ConfigurationParameter config{this};
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_GENERICPARAMETERS_H_
