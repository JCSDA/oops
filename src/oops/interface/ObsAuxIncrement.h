/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSAUXINCREMENT_H_
#define OOPS_INTERFACE_OBSAUXINCREMENT_H_

#include <iostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/ObsAuxControl.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsAuxIncrement : public util::Printable,
                        private util::ObjectCounter<ObsAuxIncrement<MODEL> > {
  typedef typename MODEL::ObsAuxIncrement     ObsAuxIncrement_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;

 public:
  static const std::string classname() {return "oops::ObsAuxIncrement";}

/// Constructor, destructor
  explicit ObsAuxIncrement(const eckit::Configuration &);
  ObsAuxIncrement(const ObsAuxIncrement &, const bool copy = true);
  ObsAuxIncrement(const ObsAuxIncrement &, const eckit::Configuration &);
  ~ObsAuxIncrement();

/// Interfacing
  const ObsAuxIncrement_ & obsauxincrement() const {return *aux_;}
  ObsAuxIncrement_ & obsauxincrement() {return *aux_;}

/// Linear algebra operators
  void diff(const ObsAuxControl_ &, const ObsAuxControl_ &);
  void zero();
  ObsAuxIncrement & operator=(const ObsAuxIncrement &);
  ObsAuxIncrement & operator+=(const ObsAuxIncrement &);
  ObsAuxIncrement & operator-=(const ObsAuxIncrement &);
  ObsAuxIncrement & operator*=(const double &);
  void axpy(const double &, const ObsAuxIncrement &);
  double dot_product_with(const ObsAuxIncrement &) const;

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<ObsAuxIncrement_> aux_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsAuxControl<MODEL> & operator+=(ObsAuxControl<MODEL> & xx, const ObsAuxIncrement<MODEL> & dx) {
  Log::trace() << "operator+=(ObsAuxControl, ObsAuxIncrement) starting" << std::endl;
  util::Timer timer("oops::ObsAuxIncrement", "operator+=ObsAuxControl");
  xx.obsauxcontrol() += dx.obsauxincrement();
  Log::trace() << "operator+=(ObsAuxControl, ObsAuxIncrement) done" << std::endl;
  return xx;
}

// =============================================================================

template<typename MODEL>
ObsAuxIncrement<MODEL>::ObsAuxIncrement(const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ObsAuxIncrement<MODEL>::ObsAuxIncrement starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxIncrement");
  aux_.reset(new ObsAuxIncrement_(conf));
  Log::trace() << "ObsAuxIncrement<MODEL>::ObsAuxIncrement done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrement<MODEL>::ObsAuxIncrement(const ObsAuxIncrement & other,
                                        const bool copy) : aux_()
{
  Log::trace() << "ObsAuxIncrement<MODEL>::ObsAuxIncrement copy starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxIncrement");
  aux_.reset(new ObsAuxIncrement_(*other.aux_, copy));
  Log::trace() << "ObsAuxIncrement<MODEL>::ObsAuxIncrement copy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrement<MODEL>::ObsAuxIncrement(const ObsAuxIncrement & other,
                                        const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ObsAuxIncrement<MODEL>::ObsAuxIncrement interpolated starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxIncrement");
  aux_.reset(new ObsAuxIncrement_(*other.aux_, conf));
  Log::trace() << "ObsAuxIncrement<MODEL>::ObsAuxIncrement interpolated done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrement<MODEL>::~ObsAuxIncrement() {
  Log::trace() << "ObsAuxIncrement<MODEL>::~ObsAuxIncrement starting" << std::endl;
  util::Timer timer(classname(), "~ObsAuxIncrement");
  aux_.reset();
  Log::trace() << "ObsAuxIncrement<MODEL>::~ObsAuxIncrement done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrement<MODEL>::diff(const ObsAuxControl_ & x1, const ObsAuxControl_ & x2) {
  Log::trace() << "ObsAuxIncrement<MODEL>::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  aux_->diff(x1.obsauxcontrol(), x2.obsauxcontrol());
  Log::trace() << "ObsAuxIncrement<MODEL>::diff done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrement<MODEL>::zero() {
  Log::trace() << "ObsAuxIncrement<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  aux_->zero();
  Log::trace() << "ObsAuxIncrement<MODEL>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrement<MODEL> & ObsAuxIncrement<MODEL>::operator=(const ObsAuxIncrement & rhs) {
  Log::trace() << "ObsAuxIncrement<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *aux_ = *rhs.aux_;
  Log::trace() << "ObsAuxIncrement<MODEL>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrement<MODEL> & ObsAuxIncrement<MODEL>::operator+=(const ObsAuxIncrement & rhs) {
  Log::trace() << "ObsAuxIncrement<MODEL>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *aux_ += *rhs.aux_;
  Log::trace() << "ObsAuxIncrement<MODEL>::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrement<MODEL> & ObsAuxIncrement<MODEL>::operator-=(const ObsAuxIncrement & rhs) {
  Log::trace() << "ObsAuxIncrement<MODEL>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *aux_ -= *rhs.aux_;
  Log::trace() << "ObsAuxIncrement<MODEL>::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsAuxIncrement<MODEL> & ObsAuxIncrement<MODEL>::operator*=(const double & zz) {
  Log::trace() << "ObsAuxIncrement<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *aux_ *= zz;
  Log::trace() << "ObsAuxIncrement<MODEL>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrement<MODEL>::axpy(const double & zz, const ObsAuxIncrement & dx) {
  Log::trace() << "ObsAuxIncrement<MODEL>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
  aux_->axpy(zz, *dx.aux_);
  Log::trace() << "ObsAuxIncrement<MODEL>::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ObsAuxIncrement<MODEL>::dot_product_with(const ObsAuxIncrement & dx) const {
  Log::trace() << "ObsAuxIncrement<MODEL>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = aux_->dot_product_with(*dx.aux_);
  Log::trace() << "ObsAuxIncrement<MODEL>::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrement<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "ObsAuxIncrement<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(conf);
  Log::trace() << "ObsAuxIncrement<MODEL>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrement<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ObsAuxIncrement<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(conf);
  Log::trace() << "ObsAuxIncrement<MODEL>::write done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ObsAuxIncrement<MODEL>::norm() const {
  Log::trace() << "ObsAuxIncrement<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = aux_->norm();
  Log::trace() << "ObsAuxIncrement<MODEL>::norm done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ObsAuxIncrement<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxIncrement<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *aux_;
  Log::trace() << "ObsAuxIncrement<MODEL>::print done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSAUXINCREMENT_H_
