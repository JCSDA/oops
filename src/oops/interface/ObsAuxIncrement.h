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
#include <memory>
#include <string>
#include <vector>

#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsAuxIncrement : public util::Printable,
                        public util::Serializable,
                        private util::ObjectCounter<ObsAuxIncrement<OBS> > {
  typedef typename OBS::ObsAuxIncrement     ObsAuxIncrement_;
  typedef ObsAuxControl<OBS>       ObsAuxControl_;

 public:
  static const std::string classname() {return "oops::ObsAuxIncrement";}

/// Constructor, destructor
  ObsAuxIncrement(const ObsSpace<OBS> &, const eckit::Configuration &);
  /// Copies \p other if \p copy is true, otherwise creates zero increment
  /// of the same size as \p other.
  ObsAuxIncrement(const ObsAuxIncrement & other, const bool copy = true);
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

/// Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<ObsAuxIncrement_> aux_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsAuxControl<OBS> & operator+=(ObsAuxControl<OBS> & xx, const ObsAuxIncrement<OBS> & dx) {
  Log::trace() << "operator+=(ObsAuxControl, ObsAuxIncrement) starting" << std::endl;
  util::Timer timer("oops::ObsAuxIncrement", "operator+=ObsAuxControl");
  xx.obsauxcontrol() += dx.obsauxincrement();
  Log::trace() << "operator+=(ObsAuxControl, ObsAuxIncrement) done" << std::endl;
  return xx;
}

// =============================================================================

template<typename OBS>
ObsAuxIncrement<OBS>::ObsAuxIncrement(const ObsSpace<OBS> & os,
                                        const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ObsAuxIncrement<OBS>::ObsAuxIncrement starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxIncrement");
  aux_.reset(new ObsAuxIncrement_(os.obsspace(), conf));
  Log::trace() << "ObsAuxIncrement<OBS>::ObsAuxIncrement done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrement<OBS>::ObsAuxIncrement(const ObsAuxIncrement & other,
                                        const bool copy) : aux_()
{
  Log::trace() << "ObsAuxIncrement<OBS>::ObsAuxIncrement copy starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxIncrement");
  aux_.reset(new ObsAuxIncrement_(*other.aux_, copy));
  Log::trace() << "ObsAuxIncrement<OBS>::ObsAuxIncrement copy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrement<OBS>::~ObsAuxIncrement() {
  Log::trace() << "ObsAuxIncrement<OBS>::~ObsAuxIncrement starting" << std::endl;
  util::Timer timer(classname(), "~ObsAuxIncrement");
  aux_.reset();
  Log::trace() << "ObsAuxIncrement<OBS>::~ObsAuxIncrement done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrement<OBS>::diff(const ObsAuxControl_ & x1, const ObsAuxControl_ & x2) {
  Log::trace() << "ObsAuxIncrement<OBS>::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  aux_->diff(x1.obsauxcontrol(), x2.obsauxcontrol());
  Log::trace() << "ObsAuxIncrement<OBS>::diff done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrement<OBS>::zero() {
  Log::trace() << "ObsAuxIncrement<OBS>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  aux_->zero();
  Log::trace() << "ObsAuxIncrement<OBS>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrement<OBS> & ObsAuxIncrement<OBS>::operator=(const ObsAuxIncrement & rhs) {
  Log::trace() << "ObsAuxIncrement<OBS>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *aux_ = *rhs.aux_;
  Log::trace() << "ObsAuxIncrement<OBS>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrement<OBS> & ObsAuxIncrement<OBS>::operator+=(const ObsAuxIncrement & rhs) {
  Log::trace() << "ObsAuxIncrement<OBS>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *aux_ += *rhs.aux_;
  Log::trace() << "ObsAuxIncrement<OBS>::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrement<OBS> & ObsAuxIncrement<OBS>::operator-=(const ObsAuxIncrement & rhs) {
  Log::trace() << "ObsAuxIncrement<OBS>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *aux_ -= *rhs.aux_;
  Log::trace() << "ObsAuxIncrement<OBS>::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
ObsAuxIncrement<OBS> & ObsAuxIncrement<OBS>::operator*=(const double & zz) {
  Log::trace() << "ObsAuxIncrement<OBS>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *aux_ *= zz;
  Log::trace() << "ObsAuxIncrement<OBS>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrement<OBS>::axpy(const double & zz, const ObsAuxIncrement & dx) {
  Log::trace() << "ObsAuxIncrement<OBS>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
  aux_->axpy(zz, *dx.aux_);
  Log::trace() << "ObsAuxIncrement<OBS>::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
double ObsAuxIncrement<OBS>::dot_product_with(const ObsAuxIncrement & dx) const {
  Log::trace() << "ObsAuxIncrement<OBS>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = aux_->dot_product_with(*dx.aux_);
  Log::trace() << "ObsAuxIncrement<OBS>::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrement<OBS>::read(const eckit::Configuration & conf) {
  Log::trace() << "ObsAuxIncrement<OBS>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(conf);
  Log::trace() << "ObsAuxIncrement<OBS>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrement<OBS>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ObsAuxIncrement<OBS>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(conf);
  Log::trace() << "ObsAuxIncrement<OBS>::write done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
double ObsAuxIncrement<OBS>::norm() const {
  Log::trace() << "ObsAuxIncrement<OBS>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = aux_->norm();
  Log::trace() << "ObsAuxIncrement<OBS>::norm done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename OBS>
size_t ObsAuxIncrement<OBS>::serialSize() const {
  Log::trace() << "ObsAuxIncrement<OBS>::serialSize" << std::endl;
  util::Timer timer(classname(), "serialSize");
  return aux_->serialSize();
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrement<OBS>::serialize(std::vector<double> & vect) const {
  Log::trace() << "ObsAuxIncrement<OBS>::serialize starting" << std::endl;
  util::Timer timer(classname(), "serialize");
  aux_->serialize(vect);
  Log::trace() << "ObsAuxIncrement<OBS>::serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrement<OBS>::deserialize(const std::vector<double> & vect, size_t & current) {
  Log::trace() << "ObsAuxIncrement<OBS>::deserialize starting" << std::endl;
  util::Timer timer(classname(), "deserialize");
  aux_->deserialize(vect, current);
  Log::trace() << "ObsAuxIncrement<OBS>::deserialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void ObsAuxIncrement<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxIncrement<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *aux_;
  Log::trace() << "ObsAuxIncrement<OBS>::print done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSAUXINCREMENT_H_
