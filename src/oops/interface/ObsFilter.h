/*
 * (C) Crown copyright 2021, Met Office
 * (C) Copyright 2024, UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "oops/base/ObsVector.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Interface class for OBS-specific implementations of the ObsFilter interface.
template <typename OBS>
class ObsFilter : public util::Printable,
                  private util::ObjectCounter<ObsVector<OBS> > {
  typedef typename OBS::ObsFilter             ObsFilter_;
  typedef GeoVaLs<OBS>                        GeoVaLs_;
  typedef ObsDiagnostics<OBS>                 ObsDiags_;
  typedef ObsSpace<OBS>                       ObsSpace_;
  typedef ObsVector<OBS>                      ObsVector_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;

 public:
  ObsFilter(const ObsSpace_ & os, const eckit::Configuration & config,
            ObsDataPtr_<int> & flags, ObsDataPtr_<float> & obserr,
            const int iteration = 0) {
    filter_ = std::make_unique<ObsFilter_>(os.obsspace(), config,
                               flags ? flags->obsdatavectorptr() : nullptr,
                               obserr ? obserr->obsdatavectorptr() : nullptr,
                               iteration);
  }

  ~ObsFilter() {
    filter_.reset();
  }

  void preProcess() {
    filter_->preProcess();
  }

  void priorFilter(const GeoVaLs<OBS> &gv) {
    filter_->priorFilter(gv.geovals());
  }

  void postFilter(const GeoVaLs<OBS> & gv, const oops::ObsVector<OBS> &ov,
                  const oops::ObsVector<OBS> &bv, const ObsDiagnostics<OBS> &dv)  {
    filter_->postFilter(gv.geovals(), ov.obsvector(), bv.obsvector(), dv.obsdiagnostics());
  }

  /// \brief Return the list of GeoVaLs required by this filter.
  Variables requiredVars() const {
    return filter_->requiredVars();
  }

  /// \brief Return the list of observation diagnostics required by this filter.
  ObsVariables requiredHdiagnostics() const {
    return filter_->requiredHdiagnostics();
  }

 private:
  void print(std::ostream & os) const override {
    os << *filter_;
  }
  std::unique_ptr<ObsFilter_> filter_;
};

// -----------------------------------------------------------------------------

}  // namespace oops
