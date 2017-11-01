/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSVECTOR_H_
#define OOPS_INTERFACE_OBSVECTOR_H_

#include <math.h>
#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/ObservationSpace.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsVector : public util::Printable,
                  private util::ObjectCounter<ObsVector<MODEL> > {
  typedef typename MODEL::ObsVector             ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsVector";}

  explicit ObsVector(const ObservationSpace<MODEL> &);
  explicit ObsVector(const ObsVector &, const bool copy = true);
  explicit ObsVector(ObsVector_ *);
  ~ObsVector();

/// Interfacing
  ObsVector_ & obsvector() {return *data_;}
  const ObsVector_ & obsvector() const {return *data_;}

// Linear algebra
  ObsVector & operator = (const ObsVector &);
  ObsVector & operator*= (const double &);
  ObsVector & operator+= (const ObsVector &);
  ObsVector & operator-= (const ObsVector &);
  ObsVector & operator*= (const ObsVector &);
  ObsVector & operator/= (const ObsVector &);

  void zero();
  void axpy(const double &, const ObsVector &);
  void invert();
  void random();
  double dot_product_with(const ObsVector &) const;
  double rms() const;
  unsigned int size() const {return data_->size();}

// I/O
  void read(const std::string &);
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<ObsVector_> data_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL>::ObsVector(const ObservationSpace<MODEL> & os): data_() {
  Log::trace() << "ObsVector<MODEL>::ObsVector starting" << std::endl;
  util::Timer timer(classname(), "ObsVector");

  data_.reset(new ObsVector_(os.observationspace()));

  Log::trace() << "ObsVector<MODEL>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL>::ObsVector(const ObsVector & other, const bool copy): data_() {
  Log::trace() << "ObsVector<MODEL>::ObsVector starting" << std::endl;
  util::Timer timer(classname(), "ObsVector");

  data_.reset(new ObsVector_(*other.data_, copy));

  Log::trace() << "ObsVector<MODEL>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL>::ObsVector(ObsVector_ * data): data_(data) {
  Log::trace() << "ObsVector<MODEL>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL>::~ObsVector() {
  Log::trace() << "ObsVector<MODEL>::~ObsVector starting" << std::endl;
  util::Timer timer(classname(), "~ObsVector");

  data_.reset();

  Log::trace() << "ObsVector<MODEL>::~ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL> & ObsVector<MODEL>::operator=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");

  *data_ = *rhs.data_;

  Log::trace() << "ObsVector<MODEL>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL> & ObsVector<MODEL>::operator*=(const double & zz) {
  Log::trace() << "ObsVector<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");

  *data_ *= zz;

  Log::trace() << "ObsVector<MODEL>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL> & ObsVector<MODEL>::operator+=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<MODEL>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");

  *data_ += *rhs.data_;

  Log::trace() << "ObsVector<MODEL>::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL> & ObsVector<MODEL>::operator-=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<MODEL>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");

  *data_ -= *rhs.data_;

  Log::trace() << "ObsVector<MODEL>::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL> & ObsVector<MODEL>::operator*=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");

  *data_ *= *rhs.data_;

  Log::trace() << "ObsVector<MODEL>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsVector<MODEL> & ObsVector<MODEL>::operator/=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<MODEL>::operator/= starting" << std::endl;
  util::Timer timer(classname(), "operator/=");

  *data_ /= *rhs.data_;

  Log::trace() << "ObsVector<MODEL>::operator/= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsVector<MODEL>::zero() {
  Log::trace() << "ObsVector<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");

  data_->zero();

  Log::trace() << "ObsVector<MODEL>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsVector<MODEL>::axpy(const double & zz, const ObsVector & rhs) {
  Log::trace() << "ObsVector<MODEL>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");

  data_->axpy(zz, *rhs.data_);

  Log::trace() << "ObsVector<MODEL>::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsVector<MODEL>::invert() {
  Log::trace() << "ObsVector<MODEL>::invert starting" << std::endl;
  util::Timer timer(classname(), "invert");

  data_->invert();

  Log::trace() << "ObsVector<MODEL>::invert done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsVector<MODEL>::random() {
  Log::trace() << "ObsVector<MODEL>::random starting" << std::endl;
  util::Timer timer(classname(), "random");

  data_->random();

  Log::trace() << "ObsVector<MODEL>::random done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
double ObsVector<MODEL>::dot_product_with(const ObsVector & other) const {
  Log::trace() << "ObsVector<MODEL>::dot_product starting" << std::endl;
  util::Timer timer(classname(), "dot_product");

  const double zz = data_->dot_product_with(*other.data_);

  Log::trace() << "ObsVector<MODEL>::dot_product done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
double ObsVector<MODEL>::rms() const {
  Log::trace() << "ObsVector<MODEL>::rms starting" << std::endl;
  util::Timer timer(classname(), "rms");

  const double zz = data_->rms();

  Log::trace() << "ObsVector<MODEL>::rms done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsVector<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsVector<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");

  os << *data_;

  Log::trace() << "ObsVector<MODEL>::print done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsVector<MODEL>::read(const std::string & name) {
  Log::trace() << "ObsVector<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");

  data_->read(name);

  Log::trace() << "ObsVector<MODEL>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsVector<MODEL>::save(const std::string & name) const {
  Log::trace() << "ObsVector<MODEL>::save starting";
  util::Timer timer(classname(), "save");

  data_->save(name);

  Log::trace() << "ObsVector<MODEL>::save done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSVECTOR_H_
