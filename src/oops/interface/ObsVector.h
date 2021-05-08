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

#include <Eigen/Dense>
#include <math.h>
#include <memory>
#include <ostream>
#include <string>

#include "oops/interface/ObsDataVector_head.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/Logger.h"
#include "oops/util/MemoryCounter.h"
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

  /// Creates vector from \p obsspace. If \p name is specified, reads the
  /// specified \p name variable from \p obsspace. Otherwise, zero vector is created.
  explicit ObsVector(const ObsSpace<OBS> & obsspace, const std::string name = "");
  ObsVector(const ObsVector &);
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

  /// Pack observations local to this MPI task into an Eigen vector
  /// (excluding vector elements that are masked out and where \p mask != 0)
  Eigen::VectorXd  packEigen(const ObsDataVector<OBS, int> & mask) const;
  /// Number of non-masked out observations local to this MPI task
  /// (size of an Eigen vector returned by `packEigen`)
  size_t packEigenSize(const ObsDataVector<OBS, int> & mask) const;

  void zero();
  /// Set this ObsVector to ones (used in tests)
  void ones();
  void axpy(const double &, const ObsVector &);
  void invert();
  void random();
  double dot_product_with(const ObsVector &) const;
  double rms() const;
  /// Mask out elements of the vector where the passed in flags are > 0
  void mask(const ObsDataVector<OBS, int> &);
  ObsVector & operator =(const ObsDataVector<OBS, float> &);

// I/O
  void save(const std::string &) const;
  void read(const std::string &);

  /// number of non-masked out observations (across all MPI tasks)
  unsigned int nobs() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsVector_> data_;
  const eckit::mpi::Comm & commTime_;
};

// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsSpace<OBS> & os, const std::string name)
  : data_(), commTime_(os.timeComm()) {
  Log::trace() << "ObsVector<OBS>::ObsVector starting " << name << std::endl;
  util::Timer timer(classname(), "ObsVector");
  util::MemoryCounter mem(classname());

  data_.reset(new ObsVector_(os.obsspace(), name));

  Log::trace() << "ObsVector<OBS>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsVector & other): data_(), commTime_(other.commTime_) {
  Log::trace() << "ObsVector<OBS>::ObsVector starting" << std::endl;
  util::Timer timer(classname(), "ObsVector");
  util::MemoryCounter mem(classname());

  data_.reset(new ObsVector_(*other.data_));

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
void ObsVector<OBS>::ones() {
  Log::trace() << "ObsVector<OBS>::ones starting" << std::endl;
  util::Timer timer(classname(), "ones");

  data_->ones();

  Log::trace() << "ObsVector<OBS>::ones done" << std::endl;
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
void ObsVector<OBS>::mask(const ObsDataVector<OBS, int> & qc) {
  Log::trace() << "ObsVector<OBS>::mask starting" << std::endl;
  util::Timer timer(classname(), "mask");
  data_->mask(qc.obsdatavector());
  Log::trace() << "ObsVector<OBS>::mask done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS> & ObsVector<OBS>::operator=(const ObsDataVector<OBS, float> & rhs) {
  Log::trace() << "ObsVector<OBS>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *data_ = rhs.obsdatavector();
  Log::trace() << "ObsVector<OBS>::operator= done" << std::endl;
  return *this;
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
template <typename OBS>
Eigen::VectorXd  ObsVector<OBS>::packEigen(const ObsDataVector<OBS, int> & mask) const {
  Log::trace() << "ObsVector<OBS>::packEigen starting " << std::endl;
  util::Timer timer(classname(), "packEigen");

  Eigen::VectorXd vec = data_->packEigen(mask.obsdatavector());

  Log::trace() << "ObsVector<OBS>::packEigen done" << std::endl;
  return vec;
}
// -----------------------------------------------------------------------------
template <typename OBS>
size_t ObsVector<OBS>::packEigenSize(const ObsDataVector<OBS, int> & mask) const {
  Log::trace() << "ObsVector<OBS>::packEigenSize starting " << std::endl;
  util::Timer timer(classname(), "packEigenSize");

  size_t len = data_->packEigenSize(mask.obsdatavector());

  Log::trace() << "ObsVector<OBS>::packEigen done" << std::endl;
  return len;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::read(const std::string & name) {
  Log::trace() << "ObsVector<OBS>::read starting " << name << std::endl;
  util::Timer timer(classname(), "read");

  data_->read(name);

  Log::trace() << "ObsVector<OBS>::read done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSVECTOR_H_
