/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSFILTER_H_
#define OOPS_INTERFACE_OBSFILTER_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
class ObsFilter : public ObsFilterBase<OBS> {
  typedef GeoVaLs<OBS>             GeoVaLs_;
  typedef ObsDiagnostics<OBS>      ObsDiags_;
  typedef ObsSpace<OBS>            ObsSpace_;
  typedef ObsVector<OBS>           ObsVector_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;
  template <typename DATA> using ObsDataVec_ = typename OBS::template ObsDataVector<DATA>;

 public:
  static const std::string classname() {return "oops::ObsFilter";}

  ObsFilter(const ObsSpace_ &, const eckit::Configuration &,
            ObsDataPtr_<int>, ObsDataPtr_<float>);
  ~ObsFilter();

  void preProcess() const override;
  void priorFilter(const GeoVaLs_ &) const override;
  void postFilter(const ObsVector_ &, const ObsDiags_ &) const override;

  Variables requiredVars() const override;
  Variables requiredHdiagnostics() const override;

 private:
  void print(std::ostream &) const override;

  ObsSpace_ obsdb_;
  const eckit::LocalConfiguration conf_;
  std::unique_ptr<FILTER> ofilt_;
};

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
ObsFilter<OBS, FILTER>::ObsFilter(const ObsSpace_ & os,
                                    const eckit::Configuration & conf,
                                    ObsDataPtr_<int> flags, ObsDataPtr_<float> obserr)
  : obsdb_(os), conf_(conf), ofilt_()
{
  Log::trace() << "ObsFilter<OBS, FILTER>::ObsFilter Configuration starting" << std::endl;
  util::Timer timer(classname(), "ObsFilter");

  std::shared_ptr<ObsDataVec_<int> > qc;
  std::shared_ptr<ObsDataVec_<float> > oberr;
  if (flags) qc = flags->obsdatavectorptr();
  if (obserr) oberr = obserr->obsdatavectorptr();

  ofilt_.reset(new FILTER(obsdb_.obsspace(), conf, qc, oberr));
  Log::trace() << "ObsFilter<OBS, FILTER>::ObsFilter Configuration done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
ObsFilter<OBS, FILTER>::~ObsFilter() {
  Log::trace() << "ObsFilter<OBS, FILTER>::~ObsFilter starting" << std::endl;
  util::Timer timer(classname(), "~ObsFilter");
  ofilt_.reset();
  Log::trace() << "ObsFilter<OBS, FILTER>::~ObsFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
void ObsFilter<OBS, FILTER>::preProcess() const {
  Log::trace() << "ObsFilter<OBS, FILTER>:: preProcess starting" << std::endl;
  util::Timer timer(classname(), "preProcess");
  ofilt_->preProcess();
  Log::trace() << "ObsFilter<OBS, FILTER>:: preProcess done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
void ObsFilter<OBS, FILTER>::priorFilter(const GeoVaLs_ & gv) const {
  Log::trace() << "ObsFilter<OBS, FILTER>:: priorFilter starting" << std::endl;
  util::Timer timer(classname(), "priorFilter");
  ofilt_->priorFilter(gv.geovals());
  Log::trace() << "ObsFilter<OBS, FILTER>:: priorFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
void ObsFilter<OBS, FILTER>::postFilter(const ObsVector_ & ov, const ObsDiags_ & dv) const {
  Log::trace() << "ObsFilter<OBS, FILTER>::postFilter starting" << std::endl;
  util::Timer timer(classname(), "postFilter");
  ofilt_->postFilter(ov.obsvector(), dv.obsdiagnostics());
  Log::trace() << "ObsFilter<OBS, FILTER>::postFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
Variables ObsFilter<OBS, FILTER>::requiredVars() const {
  Log::trace() << "ObsFilter::requiredVars" << std::endl;
  return ofilt_->requiredVars();
}

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
Variables ObsFilter<OBS, FILTER>::requiredHdiagnostics() const {
  Log::trace() << "ObsFilter::requiredHdiagnostics" << std::endl;
  return ofilt_->requiredHdiagnostics();
}

// -----------------------------------------------------------------------------

template <typename OBS, typename FILTER>
void ObsFilter<OBS, FILTER>::print(std::ostream & os) const {
  os << "ObsFilter " << conf_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSFILTER_H_
