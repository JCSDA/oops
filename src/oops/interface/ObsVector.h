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
#include <memory>
#include <ostream>
#include <string>

#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsVector : public util::Printable,
                  private util::ObjectCounter<ObsVector<OBS> > {
  typedef typename OBS::ObsVector             ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsVector";}

  explicit ObsVector(const ObsSpace<OBS> &, const std::string name = "", const bool fail = true);
  explicit ObsVector(const ObsVector &);
  ObsVector(const ObsSpace<OBS> &, const ObsVector &);
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

// accessors for data local PE
  double operator[](const size_t ii) const {return (*data_)[ii];}

  void zero();
  void axpy(const double &, const ObsVector &);
  void invert();
  void random();
  double dot_product_with(const ObsVector &) const;
  double rms() const;

// I/O
  void save(const std::string &) const;

  unsigned int nobs() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsVector_> data_;
  const eckit::mpi::Comm & commTime_;
};

// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsSpace<OBS> & os, const std::string name,
                          const bool fail): data_(), commTime_(os.timeComm()) {
  Log::trace() << "ObsVector<OBS>::ObsVector starting " << name << std::endl;
  util::Timer timer(classname(), "ObsVector");

  data_.reset(new ObsVector_(os.obsspace(), name, fail));

  Log::trace() << "ObsVector<OBS>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsVector & other): data_(), commTime_(other.commTime_) {
  Log::trace() << "ObsVector<OBS>::ObsVector starting" << std::endl;
  util::Timer timer(classname(), "ObsVector");

  data_.reset(new ObsVector_(*other.data_));

  Log::trace() << "ObsVector<OBS>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsSpace<OBS> & os, const ObsVector & other)
  : data_(), commTime_(os.timeComm())
{
  Log::trace() << "ObsVector<OBS>::ObsVector starting" << std::endl;
  util::Timer timer(classname(), "ObsVector");

  data_.reset(new ObsVector_(os.obsspace(), other.obsvector()));

  Log::trace() << "ObsVector<OBS>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::~ObsVector() {
  Log::trace() << "ObsVector<OBS>::~ObsVector starting" << std::endl;
  util::Timer timer(classname(), "~ObsVector");

  data_.reset();

  Log::trace() << "ObsVector<OBS>::~ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS> & ObsVector<OBS>::operator=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<OBS>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");

  *data_ = *rhs.data_;

  Log::trace() << "ObsVector<OBS>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS> & ObsVector<OBS>::operator*=(const double & zz) {
  Log::trace() << "ObsVector<OBS>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");

  *data_ *= zz;

  Log::trace() << "ObsVector<OBS>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS> & ObsVector<OBS>::operator+=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<OBS>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");

  *data_ += *rhs.data_;

  Log::trace() << "ObsVector<OBS>::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS> & ObsVector<OBS>::operator-=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<OBS>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");

  *data_ -= *rhs.data_;

  Log::trace() << "ObsVector<OBS>::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS> & ObsVector<OBS>::operator*=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<OBS>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");

  *data_ *= *rhs.data_;

  Log::trace() << "ObsVector<OBS>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS> & ObsVector<OBS>::operator/=(const ObsVector & rhs) {
  Log::trace() << "ObsVector<OBS>::operator/= starting" << std::endl;
  util::Timer timer(classname(), "operator/=");

  *data_ /= *rhs.data_;

  Log::trace() << "ObsVector<OBS>::operator/= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::zero() {
  Log::trace() << "ObsVector<OBS>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");

  data_->zero();

  Log::trace() << "ObsVector<OBS>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::axpy(const double & zz, const ObsVector & rhs) {
  Log::trace() << "ObsVector<OBS>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");

  data_->axpy(zz, *rhs.data_);

  Log::trace() << "ObsVector<OBS>::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::invert() {
  Log::trace() << "ObsVector<OBS>::invert starting" << std::endl;
  util::Timer timer(classname(), "invert");

  data_->invert();

  Log::trace() << "ObsVector<OBS>::invert done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::random() {
  Log::trace() << "ObsVector<OBS>::random starting" << std::endl;
  util::Timer timer(classname(), "random");

  data_->random();

  Log::trace() << "ObsVector<OBS>::random done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
double ObsVector<OBS>::dot_product_with(const ObsVector & other) const {
  Log::trace() << "ObsVector<OBS>::dot_product starting" << std::endl;
  util::Timer timer(classname(), "dot_product");

  double zz = data_->dot_product_with(*other.data_);
  commTime_.allReduceInPlace(zz, eckit::mpi::Operation::SUM);

  Log::trace() << "ObsVector<OBS>::dot_product done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template <typename OBS>
double ObsVector<OBS>::rms() const {
  Log::trace() << "ObsVector<OBS>::rms starting" << std::endl;
  util::Timer timer(classname(), "rms");

  double zz = 0.0;
  size_t ntot = this->nobs();
  if (ntot > 0) {
    zz = data_->rms();
    double zzz = zz * zz * static_cast<double>(data_->nobs());
    commTime_.allReduceInPlace(zzz, eckit::mpi::Operation::SUM);
    zzz /= static_cast<double>(ntot);
    zz = std::sqrt(zzz);
  }

  Log::trace() << "ObsVector<OBS>::rms done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template <typename OBS>
unsigned int ObsVector<OBS>::nobs() const {
  int nobs = data_->nobs();
  commTime_.allReduceInPlace(nobs, eckit::mpi::Operation::SUM);
  return nobs;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsVector<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  if (commTime_.size() > 1) {
    gatherPrint(os, *data_, commTime_);
  } else {
    os << *data_;
  }
  Log::trace() << "ObsVector<OBS>::print done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::save(const std::string & name) const {
  Log::trace() << "ObsVector<OBS>::save starting " << name << std::endl;
  util::Timer timer(classname(), "save");

  data_->save(name);

  Log::trace() << "ObsVector<OBS>::save done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSVECTOR_H_
