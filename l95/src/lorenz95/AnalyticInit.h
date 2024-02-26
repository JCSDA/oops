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

#include "lorenz95/L95Traits.h"

namespace lorenz95 {
  class LocsL95;
  class GomL95;

/// AnalyticInit class fills GeoVaLs with analytic formulae
class AnalyticInit: public oops::interface::AnalyticInitBase<L95ObsTraits> {
 public:
  explicit AnalyticInit(const eckit::Configuration &);
  void fillGeoVaLs(const LocsL95 &, GomL95 &) const override;

 private:
  eckit::LocalConfiguration config_;
};

}  // namespace lorenz95

#endif  // LORENZ95_ANALYTICINIT_H_
