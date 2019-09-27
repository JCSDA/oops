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

#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
class ObsFilter : public ObsFilterBase<MODEL> {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObsDiagnostics<MODEL>      ObsDiags_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsDataVector<MODEL, DATA> >;
  template <typename DATA> using ObsDataVec_ = typename MODEL::template ObsDataVector<DATA>;

 public:
  static const std::string classname() {return "oops::ObsFilter";}

  ObsFilter(const ObsSpace_ &, const eckit::Configuration &,
            ObsDataPtr_<int>, ObsDataPtr_<float>);
  ~ObsFilter();

  void preProcess() const override;
  void priorFilter(const GeoVaLs_ &) const override;
  void postFilter(const ObsVector_ &, const ObsDiags_ &) const override;

  Variables requiredGeoVaLs() const override;
  Variables requiredHdiagnostics() const override;

 private:
  void print(std::ostream &) const override;
  const eckit::LocalConfiguration conf_;
  std::unique_ptr<FILTER> ofilt_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
ObsFilter<MODEL, FILTER>::ObsFilter(const ObsSpace_ & os,
                                    const eckit::Configuration & conf,
                                    ObsDataPtr_<int> flags, ObsDataPtr_<float> obserr)
  : conf_(conf), ofilt_()
{
  Log::trace() << "ObsFilter<MODEL, FILTER>::ObsFilter Configuration starting" << std::endl;
  util::Timer timer(classname(), "ObsFilter");

  boost::shared_ptr<ObsDataVec_<int> > qc;
  boost::shared_ptr<ObsDataVec_<float> > oberr;
  if (flags) qc = flags->obsdatavectorptr();
  if (obserr) oberr = obserr->obsdatavectorptr();

  ofilt_.reset(new FILTER(os.observationspace(), conf, qc, oberr));
  Log::trace() << "ObsFilter<MODEL, FILTER>::ObsFilter Configuration done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
ObsFilter<MODEL, FILTER>::~ObsFilter() {
  Log::trace() << "ObsFilter<MODEL, FILTER>::~ObsFilter starting" << std::endl;
  util::Timer timer(classname(), "~ObsFilter");
  ofilt_.reset();
  Log::trace() << "ObsFilter<MODEL, FILTER>::~ObsFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
void ObsFilter<MODEL, FILTER>::preProcess() const {
  Log::trace() << "ObsFilter<MODEL, FILTER>:: preProcess starting" << std::endl;
  util::Timer timer(classname(), "preProcess");
  ofilt_->preProcess();
  Log::trace() << "ObsFilter<MODEL, FILTER>:: preProcess done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
void ObsFilter<MODEL, FILTER>::priorFilter(const GeoVaLs_ & gv) const {
  Log::trace() << "ObsFilter<MODEL, FILTER>:: priorFilter starting" << std::endl;
  util::Timer timer(classname(), "priorFilter");
  ofilt_->priorFilter(gv.geovals());
  Log::trace() << "ObsFilter<MODEL, FILTER>:: priorFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
void ObsFilter<MODEL, FILTER>::postFilter(const ObsVector_ & ov, const ObsDiags_ & dv) const {
  Log::trace() << "ObsFilter<MODEL, FILTER>::postFilter starting" << std::endl;
  util::Timer timer(classname(), "postFilter");
  ofilt_->postFilter(ov.obsvector(), dv.obsdiagnostics());
  Log::trace() << "ObsFilter<MODEL, FILTER>::postFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
Variables ObsFilter<MODEL, FILTER>::requiredGeoVaLs() const {
  Log::trace() << "ObsFilter::requiredGeoVaLs" << std::endl;
  return ofilt_->requiredGeoVaLs();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
Variables ObsFilter<MODEL, FILTER>::requiredHdiagnostics() const {
  Log::trace() << "ObsFilter::requiredHdiagnostics" << std::endl;
  return ofilt_->requiredHdiagnostics();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FILTER>
void ObsFilter<MODEL, FILTER>::print(std::ostream & os) const {
  os << "ObsFilter " << conf_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSFILTER_H_
