/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_GETVALUES_H_
#define OOPS_INTERFACE_GETVALUES_H_

#include <memory>
#include <string>


#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

/// \brief Gets values from model State to observation locations (fills GeoVaLs)
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
class GetValues : public util::Printable,
                  private util::ObjectCounter<GetValues<MODEL, OBS> > {
  typedef typename MODEL::GetValues  GetValues_;
  typedef Geometry<MODEL>            Geometry_;
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef Locations<OBS>             Locations_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::GetValues";}

/// Constructor, destructor
  GetValues(const Geometry_ &, const Locations_ &, const eckit::Configuration &);
  ~GetValues();

/// Interfacing
  GetValues_ & getvalues() {return *getvalues_;}
  const GetValues_ & getvalues() const {return *getvalues_;}

/// Get state values at observation locations
  void fillGeoVaLs(const State_ &, const util::DateTime &, const util::DateTime &,
                   GeoVaLs_ &) const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<GetValues_> getvalues_;
};

// =============================================================================

template<typename MODEL, typename OBS>
GetValues<MODEL, OBS>::GetValues(const Geometry_ & resol, const Locations_ & locs,
                                 const eckit::Configuration & conf)
    :getvalues_()
{
  Log::trace() << "GetValues<MODEL, OBS>::GetValues starting" << std::endl;
  util::Timer timer(classname(), "GetValues");
  getvalues_.reset(new GetValues_(resol.geometry(), locs.locations(), conf));
  Log::trace() << "GetValues<MODEL, OBS>::GetValues done" << std::endl;
}


// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
GetValues<MODEL, OBS>::~GetValues() {
  Log::trace() << "GetValues<MODEL, OBS>::~GetValues starting" << std::endl;
  util::Timer timer(classname(), "~GetValues");
  getvalues_.reset();
  Log::trace() << "GetValues<MODEL, OBS>::~GetValues done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::fillGeoVaLs(const State_ & state, const util::DateTime & t1,
                                   const util::DateTime & t2, GeoVaLs_ & gvals) const {
  Log::trace() << "GetValues<MODEL, OBS>::fillGeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "fillGeoVaLs");
  getvalues_->fillGeoVaLs(state.state(), t1, t2, gvals.geovals());
  Log::trace() << "GetValues<MODEL, OBS>::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::print(std::ostream & os) const {
  Log::trace() << "GetValues<MODEL, OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *getvalues_;
  Log::trace() << "GetValues<MODEL, OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_GETVALUES_H_
