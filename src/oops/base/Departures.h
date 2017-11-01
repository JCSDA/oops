/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_DEPARTURES_H_
#define OOPS_BASE_DEPARTURES_H_

#include <iostream>
#include <string>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "util/Logger.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsErrorBase.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/Printable.h"
#include "util/dot_product.h"

namespace oops {

template<typename MODEL> class Departures;

/// Difference between two observation vectors.
/*!
 * A departure is the difference between two observations.
 * The archetypal example is \f$ \mathbf{y} - {\cal H}(\mathbf{x}) \f$.
 *
 * Keeping an observation space vector here is necessary for the implementation
 * of generic observation error covariance matrices.
 */

// -----------------------------------------------------------------------------
template <typename MODEL>
class Departures : public util::Printable,
                   public GeneralizedDepartures {
  typedef ModelAtLocations<MODEL>    GOM_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncr_;
  typedef ObsErrorBase<MODEL>        ObsErrorBase_;
  typedef LinearObsOperator<MODEL>   LinearObsOperator_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
// Constructors and destructor
  explicit Departures(const ObsSpace_ &);
  explicit Departures(ObsVector_ *);
  explicit Departures(const Departures &);
  ~Departures();

// Linear algebra operators
  Departures & operator=(const Departures &);
  Departures & operator+=(const Departures &);
  Departures & operator-=(const Departures &);
  Departures & operator*=(const double &);
  Departures & operator*=(const Departures &);
  Departures & operator/=(const Departures &);
  void zero();
  void invert();
  void axpy(const double &, const Departures &);
  double dot_product_with(const Departures &) const;

/// Compute observations equivalents (TL)
  void runObsOperatorTL(const LinearObsOperator_ &, const GOM_ &, const ObsAuxIncr_ &);

/// Compute observations equivalents (AD)
  void runObsOperatorAD(const LinearObsOperator_ &, GOM_ &, ObsAuxIncr_ &) const;

/// Get departue values
  const ObsVector_ & depvalues() const {return *dep_;}

/// Save departures values
  void save(const std::string &) const;

/// Number of obs (for info only)
  unsigned int numberOfObs() const {return dep_->size();}

/// Double despatch for obs error covariance
  void helpCovarRandomize(const ObsErrorBase_ & R) {
    R.randomize(*dep_);
  }

 private:
  void print(std::ostream &) const;

  boost::scoped_ptr<ObsVector_> dep_;
};

// =============================================================================

template<typename MODEL>
Departures<MODEL>::Departures(const ObsSpace_ & obsgeom): dep_(new ObsVector_(obsgeom))
{
  Log::trace() << "Departures created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL>::Departures(ObsVector_ * val): dep_(val)
{
  Log::trace() << "Departures created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL>::Departures(const Departures & other): dep_(new ObsVector_(*other.dep_))
{
  Log::trace() << "Departures copy-created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL>::~Departures() {
  Log::trace() << "Departures destructed" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::runObsOperatorTL(const LinearObsOperator_ & hop, const GOM_ & gom,
                                         const ObsAuxIncr_ & octl) {
  hop.obsEquivTL(gom, *dep_, octl);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::runObsOperatorAD(const LinearObsOperator_ & hop, GOM_ & gom,
                                         ObsAuxIncr_ & octl) const {
  hop.obsEquivAD(gom, *dep_, octl);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator=(const Departures & rhs) {
  *dep_ = *rhs.dep_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator+=(const Departures & rhs) {
  *dep_ += *rhs.dep_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator-=(const Departures & rhs) {
  *dep_ -= *rhs.dep_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator*=(const double & zz) {
  *dep_ *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator*=(const Departures & rhs) {
  *dep_ *= *rhs.dep_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator/=(const Departures & rhs) {
  *dep_ /= *rhs.dep_;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::zero() {
  dep_->zero();
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::invert() {
  dep_->invert();
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::axpy(const double & zz, const Departures & rhs) {
  dep_->axpy(zz, *rhs.dep_);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double Departures<MODEL>::dot_product_with(const Departures & other) const {
  double zz = dot_product(*dep_, *other.dep_);
  return zz;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Departures<MODEL>::save(const std::string & name) const {
  dep_->save(name);
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Departures<MODEL>::print(std::ostream & os) const {
  os << *dep_;
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_DEPARTURES_H_
