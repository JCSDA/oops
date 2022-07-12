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
#include <utility>

#include "oops/interface/ObsDataVector_head.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

namespace interface {

// -----------------------------------------------------------------------------
/// \brief Holds observation vector (e.g. vector of observation values, or of computed H(x))
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

  /// Wraps an existing ObsVector_.
  /// This wrapping constructor doesn't need to be implemented in an ObsVector implementation.
  /// \param obsvector  The vector to wrap.
  explicit ObsVector(std::unique_ptr<ObsVector_> obsvector);
  /// Copy constructor
  ObsVector(const ObsVector &);
  /// Destructor (defined explicitly for timing and tracing)
  ~ObsVector();

  /// Accessor
  ObsVector_ & obsvector() {return *data_;}
  /// Const accessor
  const ObsVector_ & obsvector() const {return *data_;}

  /// Linear algebra operators
  ObsVector & operator = (const ObsVector &);
  ObsVector & operator*= (const double &);
  ObsVector & operator+= (const ObsVector &);
  ObsVector & operator-= (const ObsVector &);
  ObsVector & operator*= (const ObsVector &);
  ObsVector & operator/= (const ObsVector &);
  /// Add \p zz * \p rhs to the ObsVector.
  void axpy(const double & zz, const ObsVector & rhs);

  /// Pack observations local to this MPI task into an Eigen vector
  /// (excluding vector elements that are masked out: where \p mask is a missing value)
  Eigen::VectorXd packEigen(const ObsVector & mask) const;
  /// Number of non-masked out observations local to this MPI task
  /// (size of an Eigen vector returned by `packEigen`)
  size_t packEigenSize(const ObsVector & mask) const;

  /// Zero out this ObsVector
  void zero();
  /// Set this ObsVector to ones (used in tests)
  void ones();
  /// Set each value in this ObsVector to its inverse
  void invert();
  /// Set each value in this ObsVector to a random value
  void random();

  /// Return the dot product between this ObsVector and another one \p other
  double dot_product_with(const ObsVector & other) const;
  /// Return this ObsVector rms
  double rms() const;
  /// Mask out elements of the vector where \p mask is a missing value
  void mask(const ObsVector & mask);
  /// Assignment operator from \p rhs ObsDataVector<OBS, float>
  ObsVector & operator =(const ObsDataVector<OBS, float> & rhs);

  /// Save this ObsVector as group \p name in the ObsSpace
  void save(const std::string &) const;
  /// Fill ObsVector with data with group \p name from the associated ObsSpace
  void read(const std::string &);

  /// Number of non-masked out observations (across all MPI tasks)
  unsigned int nobs() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsVector_> data_;
};

// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsSpace<OBS> & os, const std::string name) : data_() {
  Log::trace() << "ObsVector<OBS>::ObsVector starting " << name << std::endl;
  util::Timer timer(classname(), "ObsVector");
  data_.reset(new ObsVector_(os.obsspace(), name));
  this->setObjectSize(data_->size() * sizeof(double));
  Log::trace() << "ObsVector<OBS>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::ObsVector(std::unique_ptr<ObsVector_> obsvector)
  : data_(std::move(obsvector)) {
  Log::trace() << "ObsVector<OBS>::ObsVector starting " << std::endl;
  util::Timer timer(classname(), "ObsVector");
  this->setObjectSize(data_->size() * sizeof(double));
  Log::trace() << "ObsVector<OBS>::ObsVector done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename OBS>
ObsVector<OBS>::ObsVector(const ObsVector & other): data_() {
  Log::trace() << "ObsVector<OBS>::ObsVector starting" << std::endl;
  util::Timer timer(classname(), "ObsVector");
  data_.reset(new ObsVector_(*other.data_));
  this->setObjectSize(data_->size() * sizeof(double));
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

  Log::trace() << "ObsVector<OBS>::dot_product done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::mask(const ObsVector & mask) {
  Log::trace() << "ObsVector<OBS>::mask(ObsVector) starting" << std::endl;
  util::Timer timer(classname(), "mask(ObsVector)");
  data_->mask(mask.obsvector());
  Log::trace() << "ObsVector<OBS>::mask(ObsVector) done" << std::endl;
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

  double zz = data_->rms();

  Log::trace() << "ObsVector<OBS>::rms done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template <typename OBS>
unsigned int ObsVector<OBS>::nobs() const {
  int nobs = data_->nobs();
  return nobs;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void ObsVector<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsVector<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *data_;
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
Eigen::VectorXd  ObsVector<OBS>::packEigen(const ObsVector & mask) const {
  Log::trace() << "ObsVector<OBS>::packEigen starting " << std::endl;
  util::Timer timer(classname(), "packEigen");

  Eigen::VectorXd vec = data_->packEigen(mask.obsvector());

  Log::trace() << "ObsVector<OBS>::packEigen done" << std::endl;
  return vec;
}
// -----------------------------------------------------------------------------
template <typename OBS>
size_t ObsVector<OBS>::packEigenSize(const ObsVector & mask) const {
  Log::trace() << "ObsVector<OBS>::packEigenSize starting " << std::endl;
  util::Timer timer(classname(), "packEigenSize");

  size_t len = data_->packEigenSize(mask.obsvector());

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

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSVECTOR_H_
