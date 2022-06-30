/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_GEOVALSWRITER_H_
#define OOPS_GENERIC_GEOVALSWRITER_H_

#include <memory>

#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/generic/ObsFilterBase.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {

/// \brief Parameters controlling GeoVaLsWriter
template <typename OBS>
class GeoVaLsWriterParameters : public ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(GeoVaLsWriterParameters, ObsFilterParametersBase)
  // GeoVaLsWriter uses GeoVaLs Parameters (used in I/O)
  typedef typename GeoVaLs<OBS>::Parameters_   GeoVaLsParameters_;
 public:
  GeoVaLsParameters_ geovals{this};
};

// -----------------------------------------------------------------------------

template <typename OBS>
class GeoVaLsWriter : public ObsFilterBase<OBS> {
  typedef GeoVaLs<OBS>             GeoVaLs_;
  typedef ObsDiagnostics<OBS>      ObsDiags_;
  typedef ObsSpace<OBS>            ObsSpace_;
  typedef ObsVector<OBS>           ObsVector_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;

 public:
  typedef GeoVaLsWriterParameters<OBS> Parameters_;

  GeoVaLsWriter(const ObsSpace_ &, const Parameters_ & params,
                ObsDataPtr_<int>, ObsDataPtr_<float>): params_(params), novars_() {}
  ~GeoVaLsWriter() = default;

  void preProcess() override {}

  void priorFilter(const GeoVaLs_ & gv) override {
    const double zz = sqrt(dot_product(gv, gv));
    Log::info() << "GeoVaLsWriter norm = " << zz << std::endl;
    gv.write(params_.geovals);
  }

  void postFilter(const GeoVaLs_ & gv,
                  const ObsVector_ &,
                  const ObsVector_ &,
                  const ObsDiags_ &) override {
    this->priorFilter(gv);
  }
  void checkFilterData(const FilterStage filterStage) override {}

  Variables requiredVars() const override {return novars_;};
  Variables requiredHdiagnostics() const override {return novars_;};

 private:
  const Parameters_ params_;
  const Variables novars_;  // could be used to determine what needs saving

  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

template <typename OBS>
void GeoVaLsWriter<OBS>::print(std::ostream & os) const {
  os << "Filter outputting GeoVaLs" << params_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_GEOVALSWRITER_H_
