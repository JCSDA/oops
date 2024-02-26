/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_ANALYTICINIT_H_
#define QG_MODEL_ANALYTICINIT_H_

#include "eckit/config/LocalConfiguration.h"

#include "oops/interface/AnalyticInitBase.h"

#include "oops/qg/QgTraits.h"

namespace qg {
  class LocationsQG;
  class GomQG;

/// AnalyticInit class fills GeoVaLs with analytic formulae
/// Options: baroclinic instability and large vortices
class AnalyticInit : public oops::interface::AnalyticInitBase<QgObsTraits> {
 public:
  explicit AnalyticInit(const eckit::Configuration &);
  void fillGeoVaLs(const LocationsQG &, GomQG &) const override;

 private:
  const eckit::LocalConfiguration config_;
};

}  // namespace qg

#endif  // QG_MODEL_ANALYTICINIT_H_
