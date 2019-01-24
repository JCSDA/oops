/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_BASE_GEOVALSWRITER_H_
#define OOPS_BASE_GEOVALSWRITER_H_

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class GeoVaLsWriter : public ObsFilterBase<MODEL> {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  GeoVaLsWriter(const ObsSpace_ &, const eckit::Configuration & conf): conf_(conf), geovars_() {}
  ~GeoVaLsWriter() {}

  void priorFilter(const GeoVaLs_ & gv) const override {
    const double zz = sqrt(dot_product(gv, gv));
    Log::debug() << "GeoVaLsWriter norm = " << zz << std::endl;
    gv.write(conf_);
  }

  void postFilter(const ObsVector_ &) const override {}

  const Variables & requiredGeoVaLs() const override {return geovars_;};

 private:
  const eckit::LocalConfiguration conf_;
  const Variables geovars_;  // could be used to determine what needs saving
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void GeoVaLsWriter<MODEL>::print(std::ostream & os) const {
  os << "GeoVaLsWriter: " << conf_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GEOVALSWRITER_H_
