/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVATIONS_H_
#define OOPS_BASE_OBSERVATIONS_H_

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "oops/base/Departures.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "util/DateTime.h"
#include "util/Printable.h"

namespace oops {

template<typename MODEL> class Observations;

/// Observations Class.
/*!
 *  Contains observed values or their model equivalents and the associated
 *  observation operator.
 */

// -----------------------------------------------------------------------------
template <typename MODEL> class Observations : public util::Printable {
  typedef Departures<MODEL>          Departures_;
  typedef ModelAtLocations<MODEL>    GOM_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  explicit Observations(const ObsSpace_ &);
  explicit Observations(const Observations &);
  ~Observations();
  Observations & operator=(const Observations &);

/// Interactions with Departures
  ObsVector_ * operator-(const Observations & other) const;
  Observations & operator+=(const Departures_ &);

/// Compute observations equivalents
  void runObsOperator(const ObsOperator_ &, const GOM_ &, const ObsAuxCtrl_ &);

/// Get observations values
  const ObsVector_ & obsvalues() const {return obs_;}

/// Save observations values
  void save(const std::string &) const;
  void read(const eckit::Configuration &);

 private:
  void print(std::ostream &) const;

  ObsVector_ obs_;
};

// =============================================================================

template <typename MODEL>
Observations<MODEL>::Observations(const ObsSpace_ & obsdb): obs_(obsdb)
{
  Log::trace() << "Observations created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::Observations(const Observations & other): obs_(other.obs_)
{
  Log::trace() << "Observations copy-created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::~Observations() {
  Log::trace() << "Observations destructed" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator=(const Observations & rhs) {
  obs_ = rhs.obs_;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL> * Observations<MODEL>::operator-(const Observations & other) const {
  ObsVector_ * ovec = new ObsVector_(obs_, true);
  *ovec -= other.obs_;
  return ovec;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator+=(const Departures_ & dy) {
  obs_ += dy.depvalues();
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::runObsOperator(const ObsOperator_ & hop, const GOM_ & gom,
                                         const ObsAuxCtrl_ & ybias) {
  hop.obsEquiv(gom, obs_, ybias);
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::save(const std::string & name) const {
  obs_.save(name);
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::read(const eckit::Configuration & config) {
  const std::string name = config.getString("ObsData.obsvalue");
  obs_.read(name);
  Log::trace() << "Observations:Observations have been read" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::print(std::ostream & os) const {
  os << obs_;
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_OBSERVATIONS_H_
