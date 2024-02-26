/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_ANALYTICINIT_H_
#define OOPS_BASE_ANALYTICINIT_H_

#include <memory>
#include <string>

#include "oops/generic/AnalyticInitBase.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Timer.h"

namespace oops {

/// \brief Initializes GeoVaLs with analytic formula
template <typename OBS>
class AnalyticInit : private util::ObjectCounter<AnalyticInit<OBS> > {
  typedef AnalyticInitBase<OBS> AnalyticInit_;
  typedef GeoVaLs<OBS>          GeoVaLs_;
  typedef SampledLocations<OBS> SampledLocations_;

 public:
  static const std::string classname() {return "oops::AnalyticInit";}

  /// constructor (parameters)
  explicit AnalyticInit(const eckit::Configuration &);

  /// destructor and copy/move constructors/assignments
  ~AnalyticInit();
  AnalyticInit(const AnalyticInit &) = delete;
  AnalyticInit(AnalyticInit &&) = delete;
  AnalyticInit& operator=(const AnalyticInit &) = delete;
  AnalyticInit& operator=(AnalyticInit &&) = delete;

  /// Fill GeoVaLs with values computed by analytic function
  void fillGeoVaLs(const SampledLocations_ &, GeoVaLs_ &) const;

 private:
  std::unique_ptr<AnalyticInit_> analytic_;
};

// -----------------------------------------------------------------------------

template<typename OBS>
AnalyticInit<OBS>::AnalyticInit(const eckit::Configuration & config) {
  Log::trace() << "AnalyticInit<OBS>::AnalyticInit starting" << std::endl;
  util::Timer timer(classname(), "AnalyticInit");
  analytic_ = AnalyticInitFactory<OBS>::create(config);
  Log::trace() << "AnalyticInit<OBS>::AnalyticInit done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
AnalyticInit<OBS>::~AnalyticInit() {
  Log::trace() << "AnalyticInit<OBS>::~AnalyticInit starting" << std::endl;
  util::Timer timer(classname(), "~AnalyticInit");
  analytic_.reset();
  Log::trace() << "AnalyticInit<OBS>::~AnalyticInit done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void AnalyticInit<OBS>::fillGeoVaLs(const SampledLocations_ & locs,
                                    GeoVaLs_ & gvals) const {
  Log::trace() << "AnalyticInit<OBS>::fillGeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "fillGeoVaLs");
  analytic_->fillGeoVaLs(locs, gvals);
  Log::trace() << "AnalyticInit<OBS>::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_ANALYTICINIT_H_
