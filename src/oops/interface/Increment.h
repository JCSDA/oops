/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_INCREMENT_H_
#define OOPS_INTERFACE_INCREMENT_H_

#include <string>

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/PostProcessorTL.h"
#include "oops/base/PostProcessorAD.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/interface/State.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace oops {

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in an Increment.
 */

// -----------------------------------------------------------------------------

template <typename MODEL>
class Increment : public oops::GeneralizedDepartures,
                  public util::Printable,
                  private util::ObjectCounter<Increment<MODEL> > {
  typedef typename MODEL::Increment  Increment_;
  typedef Geometry<MODEL>            Geometry_;
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef Locations<MODEL>           Locations_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::Increment";}

/// Constructor, destructor
  Increment(const Geometry_ &, const Variables &, const util::DateTime &);
  Increment(const Geometry_ &, const Increment &);
  Increment(const Increment &, const bool copy = true);
  virtual ~Increment();

/// Interfacing
  Increment_ & increment() {return *increment_;}
  const Increment_ & increment() const {return *increment_;}

/// Interpolate to observation location
  void interpolateTL(const Locations_ &, const Variables &, GeoVaLs_ &) const;
  void interpolateAD(const Locations_ &, const Variables &, const GeoVaLs_ &);

/// Interactions with State
  void diff(const State_ &, const State_ &);

/// Time
  const util::DateTime validTime() const {return increment_->validTime();}
  void updateTime(const util::Duration & dt) {increment_->updateTime(dt);}

/// Linear algebra operators
  void zero();
  void zero(const util::DateTime &);
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  void axpy(const double &, const Increment &, const bool check = true);
  double dot_product_with(const Increment &) const;
  void schur_product_with(const Increment &);
  void random();
  void accumul(const double &, const State_ &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

/// Get geometry
  Geometry_ geometry() const;

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<Increment_> increment_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
State<MODEL> & operator+=(State<MODEL> & xx, const Increment<MODEL> & dx) {
  Log::trace() << "operator+=(State, Increment) starting" << std::endl;
  util::Timer timer("oops::Increment", "operator+=(State, Increment)");
  xx.state() += dx.increment();
  Log::trace() << "operator+=(State, Increment) done" << std::endl;
  return xx;
}

// =============================================================================
/// Constructor, destructor
// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & resol, const Variables & vars,
                            const util::DateTime & time) : increment_()
{
  Log::trace() << "Increment<MODEL>::Increment starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(resol.geometry(), vars, time));
  Log::trace() << "Increment<MODEL>::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & resol, const Increment & other)
  : increment_()
{
  Log::trace() << "Increment<MODEL>::Increment starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(resol.geometry(), *other.increment_));
  Log::trace() << "Increment<MODEL>::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Increment & other, const bool copy) : increment_()
{
  Log::trace() << "Increment<MODEL>::Increment copy starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(*other.increment_, copy));
  Log::trace() << "Increment<MODEL>::Increment copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::~Increment() {
  Log::trace() << "Increment<MODEL>::~Increment starting" << std::endl;
  util::Timer timer(classname(), "~Increment");
  increment_.reset();
  Log::trace() << "Increment<MODEL>::~Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::interpolateTL(const Locations_ & loc, const Variables & vars,
                                     GeoVaLs_ & gvals) const {
  Log::trace() << "Increment<MODEL>::interpolateTL starting" << std::endl;
  util::Timer timer(classname(), "interpolateTL");
  increment_->interpolateTL(loc.locations(), vars, gvals.geovals());
  Log::trace() << "Increment<MODEL>::interpolateTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::interpolateAD(const Locations_ & loc, const Variables & vars,
                                     const GeoVaLs_ & gvals) {
  Log::trace() << "Increment<MODEL>::interpolateAD starting" << std::endl;
  util::Timer timer(classname(), "interpolateAD");
  increment_->interpolateAD(loc.locations(), vars, gvals.geovals());
  Log::trace() << "Increment<MODEL>::interpolateAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::diff(const State_ & x1, const State_ & x2) {
  Log::trace() << "Increment<MODEL>::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  increment_->diff(x1.state(), x2.state());
  Log::trace() << "Increment<MODEL>::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::zero() {
  Log::trace() << "Increment<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  increment_->zero();
  Log::trace() << "Increment<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::zero(const util::DateTime & tt) {
  Log::trace() << "Increment<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  increment_->zero(tt);
  Log::trace() << "Increment<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *increment_ = *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator+=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *increment_ += *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator-=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *increment_ -= *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator*=(const double & zz) {
  Log::trace() << "Increment<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *increment_ *= zz;
  Log::trace() << "Increment<MODEL>::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::axpy(const double & zz, const Increment & dx, const bool check) {
  Log::trace() << "Increment<MODEL>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
  increment_->axpy(zz, *dx.increment_, check);
  Log::trace() << "Increment<MODEL>::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::dot_product_with(const Increment & dx) const {
  Log::trace() << "Increment<MODEL>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = increment_->dot_product_with(*dx.increment_);
  Log::trace() << "Increment<MODEL>::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::schur_product_with(const Increment & dx) {
  Log::trace() << "Increment<MODEL>::schur_product_with starting" << std::endl;
  util::Timer timer(classname(), "schur_product_with");
  increment_->schur_product_with(*dx.increment_);
  Log::trace() << "Increment<MODEL>::schur_product_with done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::random() {
  Log::trace() << "Increment<MODEL>::random starting" << std::endl;
  util::Timer timer(classname(), "random");
  increment_->random();
  Log::trace() << "Increment<MODEL>::random done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::accumul(const double & zz, const State_ & xx) {
  Log::trace() << "Increment<MODEL>::accumul starting" << std::endl;
  util::Timer timer(classname(), "accumul");
  increment_->accumul(zz, xx.state());
  Log::trace() << "Increment<MODEL>::accumul done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "Increment<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  increment_->read(conf);
  Log::trace() << "Increment<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "Increment<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  increment_->write(conf);
  Log::trace() << "Increment<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::norm() const {
  Log::trace() << "Increment<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = increment_->norm();
  Log::trace() << "Increment<MODEL>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Geometry<MODEL> Increment<MODEL>::geometry() const {
  Log::trace() << "Increment<MODEL>::geometry starting" << std::endl;
  util::Timer timer(classname(), "geometry");
  Geometry<MODEL> geom(increment_->geometry());
  Log::trace() << "Increment<MODEL>::geometry done" << std::endl;
  return geom;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Increment<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *increment_;
  Log::trace() << "Increment<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_INCREMENT_H_
