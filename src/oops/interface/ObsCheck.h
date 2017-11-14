/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_INTERFACE_OBSCHECK_H_
#define OOPS_INTERFACE_OBSCHECK_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/FilterBase.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "util/dot_product.h"
#include "util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsCheck : public FilterBase<MODEL> {
  typedef typename MODEL::ObsCheck   ObsCheck_;
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsCheck";}

  ObsCheck(const ObsSpace_ &);
  ObsCheck(const eckit::Configuration & conf);
  ~ObsCheck();

  void postFilter(const GeoVaLs_ &, const ObsVector_ &, const ObsSpace_ &) const;
  void priorFilter(const GeoVaLs_ &, const ObsVector_ &, const ObsSpace_ &) const;

 private:
  const eckit::LocalConfiguration conf_;
  void print(std::ostream &) const override;
  boost::scoped_ptr<ObsCheck_> obsc_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsCheck<MODEL>::ObsCheck(const eckit::Configuration & conf) : conf_(conf) {
  Log::trace() << "ObsCheck<MODEL>::ObsCheck Configuration starting" << std::endl;
  util::Timer timer(classname(), "ObsCheck");
    obsc_.reset(new ObsCheck_(conf));
  Log::trace() << "ObsCheck<MODEL>::ObsCheck Configuration done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsCheck<MODEL>::ObsCheck(const ObsSpace_ & os) : obsc_() {
  Log::trace() << "ObsCheck<MODEL>::ObsCheck ObsSpace starting" << std::endl;
  util::Timer timer(classname(), "ObsCheck");
  Log::trace() << "ObsCheck<MODEL>::ObsCheck ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsCheck<MODEL>::~ObsCheck() {
  Log::trace() << "ObsCheck<MODEL>::~ObsCheck starting" << std::endl;
  util::Timer timer(classname(), "~ObsCheck");
    obsc_.reset();
  Log::trace() << "ObsCheck<MODEL>::~ObsCheck done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsCheck<MODEL>::postFilter(const GeoVaLs_ & gv, const ObsVector_ & ov, const ObsSpace_ & os) const {
  Log::trace() << "ObsCheck<MODEL>::postFilter starting" << std::endl;
  util::Timer timer(classname(), "ObsCheck");
    obsc_->postFilter(gv.geovals(),ov.obsvector(),os.observationspace());
  Log::trace() << "ObsCheck<MODEL>::postFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsCheck<MODEL>::priorFilter(const GeoVaLs_ & gv, const ObsVector_ & ov, const ObsSpace_ & os) const {
  Log::trace() << "ObsCheck<MODEL>:: priorFilter starting" << std::endl;
  util::Timer timer(classname(), "ObsCheck");
    obsc_->priorFilter(gv.geovals(),ov.obsvector(),os.observationspace());
  Log::trace() << "ObsCheck<MODEL>:: priorFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsCheck<MODEL>::print(std::ostream & os) const {
  os << "ObsCheck " << conf_ ;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSCHECK_H_
