/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_POSTTIMERPARAMETERS_H_
#define OOPS_BASE_POSTTIMERPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Options controlling PostTimer
class PostTimerParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(PostTimerParameters, Parameters)
 public:
  /// frequency of calling the PostProcessor (default = 0 -- call at every step)
  oops::Parameter<util::Duration> frequency{"frequency", util::Duration(0), this};
  /// constrols delta for the first call of PostProcessor (first call will happen at begin+first)
  oops::Parameter<util::Duration> first{"first", util::Duration(0), this};
  /// specifies at which times to call PostProcessor
  oops::Parameter<std::vector<util::DateTime>> times{"times", std::vector<util::DateTime>(), this};
  oops::Parameter<std::vector<util::Duration>> steps{"steps", std::vector<util::Duration>(), this};
};

}  // namespace oops

#endif  // OOPS_BASE_POSTTIMERPARAMETERS_H_
