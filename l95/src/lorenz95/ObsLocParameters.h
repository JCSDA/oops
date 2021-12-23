/*
 * (C) Copyright 2021- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSLOCPARAMETERS_H_
#define LORENZ95_OBSLOCPARAMETERS_H_

#include "oops/base/ObsLocalizationParametersBase.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace lorenz95 {

/// \brief Parameters controlling obs-space localization.
class ObsLocParameters : public oops::ObsLocalizationParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsLocParameters, oops::ObsLocalizationParametersBase)

 public:
  oops::RequiredParameter<double> lengthscale{"lengthscale",
        "Localization distance (distance where localization goes to zero)", this};
};

}  // namespace lorenz95

#endif  // LORENZ95_OBSLOCPARAMETERS_H_
