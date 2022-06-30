/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_ANALYTICINIT_H_
#define LORENZ95_ANALYTICINIT_H_

#include "eckit/config/LocalConfiguration.h"

#include "oops/interface/AnalyticInitBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

#include "lorenz95/L95Traits.h"

namespace lorenz95 {

/// Parameters for L95 Analytic init
class AnalyticInitParameters : public oops::AnalyticInitParametersBase {
  OOPS_CONCRETE_PARAMETERS(AnalyticInitParameters, AnalyticInitParametersBase)

 public:
  /// mean of the generated values
  oops::Parameter<double> mean{"mean", 0.0, this};

  /// sinus function added to the mean
  oops::OptionalParameter<double> sinus{"sinus", this};
};

// Forward declarations
class LocsL95;
class GomL95;

/// AnalyticInit class fills GeoVaLs with analytic formulae
class AnalyticInit: public oops::interface::AnalyticInitBase<L95ObsTraits> {
 public:
  typedef AnalyticInitParameters Parameters_;

  explicit AnalyticInit(const Parameters_ &);
  void fillGeoVaLs(const LocsL95 &, GomL95 &) const override;

 private:
  Parameters_ options_;
};

}  // namespace lorenz95

#endif  // LORENZ95_ANALYTICINIT_H_
