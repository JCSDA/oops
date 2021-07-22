/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_LINEARGETVALUES_H_
#define OOPS_INTERFACE_LINEARGETVALUES_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

/// \brief sets trajectory and computes TL and AD for GetValues
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
class LinearGetValues : public util::Printable,
                        private util::ObjectCounter<LinearGetValues<MODEL, OBS> > {
  typedef typename MODEL::LinearGetValues  LinearGetValues_;
  typedef Geometry<MODEL>            Geometry_;
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef Increment<MODEL>           Increment_;
  typedef Locations<OBS>             Locations_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::LinearGetValues";}

/// Constructor, destructor
  LinearGetValues(const Geometry_ &, const Locations_ &, const eckit::Configuration &);
  virtual ~LinearGetValues();

/// Interfacing
  LinearGetValues_ & lingetvalues() {return *lingetvalues_;}
  const LinearGetValues_ & lingetvalues() const {return *lingetvalues_;}

/// set trajectory for GetValues
  void setTrajectory(const State_ &, const util::DateTime &, const util::DateTime &,
                     GeoVaLs_ &);
/// compute TL of GetValues
  void fillGeoVaLsTL(const Increment_ &, const util::DateTime &, const util::DateTime &,
                     GeoVaLs_ &) const;
/// compute AD of GetValues
  void fillGeoVaLsAD(Increment_ &, const util::DateTime &, const util::DateTime &,
                     const GeoVaLs_ &) const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<LinearGetValues_> lingetvalues_;
};

// =============================================================================
/// Constructor, destructor
// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
LinearGetValues<MODEL, OBS>::LinearGetValues(const Geometry_ & resol,
                                        const Locations_ & loc,
                                        const eckit::Configuration & linearGetValuesConf)
     : lingetvalues_() {
  Log::trace() << "LinearGetValues<MODEL, OBS>::LinearGetValues starting" << std::endl;
  util::Timer timer(classname(), "LinearGetValues");
  lingetvalues_.reset(new LinearGetValues_(resol.geometry(), loc.locations(),
                                           linearGetValuesConf));
  Log::trace() << "LinearGetValues<MODEL, OBS>::LinearGetValues done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
LinearGetValues<MODEL, OBS>::~LinearGetValues() {
  Log::trace() << "LinearGetValues<MODEL, OBS>::~LinearGetValues starting" << std::endl;
  util::Timer timer(classname(), "~LinearGetValues");
  lingetvalues_.reset();
  Log::trace() << "LinearGetValues<MODEL, OBS>::~LinearGetValues done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void LinearGetValues<MODEL, OBS>::setTrajectory(const State_ & state, const util::DateTime & t1,
                                           const util::DateTime & t2, GeoVaLs_ & gvals) {
  Log::trace() << "LinearGetValues<MODEL, OBS>::setTrajectory starting" << std::endl;
  util::Timer timer(classname(), "setTrajectory");
  lingetvalues_->setTrajectory(state.state(), t1, t2, gvals.geovals());
  Log::trace() << "LinearGetValues<MODEL, OBS>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void LinearGetValues<MODEL, OBS>::fillGeoVaLsTL(const Increment_ & inc, const util::DateTime & t1,
                                           const util::DateTime & t2, GeoVaLs_ & gvals) const {
  Log::trace() << "LinearGetValues<MODEL, OBS>::fillGeoVaLsTL starting" << std::endl;
  util::Timer timer(classname(), "fillGeoVaLsTL");
  lingetvalues_->fillGeoVaLsTL(inc.increment(), t1, t2, gvals.geovals());
  Log::trace() << "LinearGetValues<MODEL, OBS>::fillGeoVaLsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void LinearGetValues<MODEL, OBS>::fillGeoVaLsAD(Increment_ & inc, const util::DateTime & t1,
                                     const util::DateTime & t2, const GeoVaLs_ & gvals) const {
  Log::trace() << "LinearGetValues<MODEL, OBS>::fillGeoVaLsAD starting" << std::endl;
  util::Timer timer(classname(), "fillGeoVaLsAD");
  lingetvalues_->fillGeoVaLsAD(inc.increment(), t1, t2, gvals.geovals());
  Log::trace() << "LinearGetValues<MODEL, OBS>::fillGeoVaLsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void LinearGetValues<MODEL, OBS>::print(std::ostream & os) const {
  Log::trace() << "LinearGetValues<MODEL, OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *lingetvalues_;
  Log::trace() << "LinearGetValues<MODEL, OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEARGETVALUES_H_
