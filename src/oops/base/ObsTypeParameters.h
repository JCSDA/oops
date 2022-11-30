/*
 * (C) Copyright 2020-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/Observer.h"
#include "oops/generic/ObsErrorBase.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/algorithms.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Options controlling the processing of observations from a single obs space.
///
/// This is an abstract base class; it can be inherited from and extended with extra parameters.
template <typename OBS>
class ObsTypeParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(ObsTypeParametersBase, Parameters)

 public:
  typedef typename ObsAuxControl<OBS>::Parameters_ ObsAuxControlParameters_;
  typedef ObsErrorParametersWrapper<OBS> ObsErrorParameters_;
  typedef typename ObsSpace<OBS>::Parameters_ ObsSpaceParameters_;

  /// Options used to configure the observation space.
  oops::RequiredParameter<ObsSpaceParameters_> obsSpace{"obs space", this};

  /// Options used to configure the observation operator, observation filters and GetValues.
  ObserverParameters<OBS> observer{this};

  /// Options used to configure the observation error covariance matrix model.
  oops::Parameter<ObsErrorParameters_> obsError{"obs error", {}, this};

  /// Options used to configure bias correction.
  oops::Parameter<ObsAuxControlParameters_> obsBias{"obs bias", {}, this};
};

// -----------------------------------------------------------------------------

/// \brief Options controlling the processing of observations from a single obs space.
template <typename OBS>
class ObsTypeParameters : public ObsTypeParametersBase<OBS> {
  OOPS_CONCRETE_PARAMETERS(ObsTypeParameters, ObsTypeParametersBase<OBS>)
};

// -----------------------------------------------------------------------------

template <typename OBS>
std::vector<typename ObsSpace<OBS>::Parameters_> obsSpaceParameters(
    const std::vector<ObsTypeParameters<OBS>> &observations) {
  return util::transformVector(
    observations,
    [](const ObsTypeParameters<OBS> & obsTypeParams) { return obsTypeParams.obsSpace.value(); });
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::vector<typename ObsAuxControl<OBS>::Parameters_> obsAuxParameters(
    const std::vector<ObsTypeParameters<OBS>> &observations) {
  return util::transformVector(
    observations,
    [](const ObsTypeParameters<OBS> & obsTypeParams) { return obsTypeParams.obsBias.value(); });
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::vector<ObsErrorParametersWrapper<OBS>> obsErrorParameters(
    const std::vector<ObsTypeParameters<OBS>> &observations) {
  return util::transformVector(
    observations,
    [](const ObsTypeParameters<OBS> & obsTypeParams) { return obsTypeParams.obsError.value(); });
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::vector<ObserverParameters<OBS>> observerParameters(
    const std::vector<ObsTypeParameters<OBS>> &observations) {
  return util::transformVector(
    observations,
    [](const ObsTypeParameters<OBS> & obsTypeParams) { return obsTypeParams.observer; });
}

// -----------------------------------------------------------------------------

}  // namespace oops

