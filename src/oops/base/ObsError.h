/*
 * (C) Crown copyright 2021, Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSERROR_H_
#define OOPS_BASE_OBSERROR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/ObsVector.h"
#include "oops/generic/ObsErrorBase.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Observation error covariance matrix of observations from a single ObsSpace.
template<typename OBS>
class ObsError : public util::Printable,
                 private util::ObjectCounter<ObsError<OBS> > {
  typedef ObsErrorBase<OBS>     ObsErrorBase_;
  typedef ObsVector<OBS>        ObsVector_;
  typedef ObsSpace<OBS>         ObsSpace_;

 public:
  static const std::string classname() {return "oops::ObsError";}

  ObsError(const ObsErrorParametersBase & params, const ObsSpace_ & os);
  ~ObsError() override;
  ObsError(const ObsError &) = delete;
  ObsError(ObsError &&) = default;
  ObsError& operator=(const ObsError &) = delete;
  ObsError& operator=(ObsError &&) = default;

  /// Multiply a Departure \p dy by \f$R\f$.
  void multiply(ObsVector_ & dy) const;
  /// Multiply a Departure \p dy by \f$R^{-1}\f$.
  void inverseMultiply(ObsVector_ & dy) const;

  /// Multiply a Departure \p dy by \f$R^{1/2}\f$.
  void sqrtMultiply(ObsVector_ & dy) const;
  /// Multiply a Departure \p dy by \f$R^{-1/2}\f$.
  void invSqrtMultiply(ObsVector_ & dy) const;

  /// Generate a random perturbation in \p dy.
  void randomize(ObsVector_ & dy) const;

  /// Save obs errors.
  void save(const std::string &) const;

  /// Return a copy of obs error std. dev. If this ObsVector_ is modified (e.g. by obs filters),
  /// it should be passed back to update() to ensure the covariance matrix stays consistent.
  ObsVector_ obserrors() const;

  /// Set the diagonal of the covariance matrix to \p stddev squared.
  void update(const ObsVector_ &stddev);

  /// Return the vector of inverse obs error variances.
  ObsVector_ inverseVariance() const;

  /// Get mean error for Jo table.
  double getRMSE() const;

 private:
  void print(std::ostream &) const override;

 private:
  std::unique_ptr<ObsErrorBase_> err_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsError<OBS>::ObsError(const ObsErrorParametersBase & params, const ObsSpace_ & os) {
  Log::trace() << "ObsError<OBS>::ObsError starting" << std::endl;

  util::Timer timer(classname(), "ObsErrors");
  size_t init = eckit::system::ResourceUsage().maxResidentSetSize();

  err_ = ObsErrorFactory<OBS>::create(params, os);

  size_t current = eckit::system::ResourceUsage().maxResidentSetSize();
  this->setObjectSize(current - init);
  Log::trace() << "ObsError<OBS>::ObsError done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsError<OBS>::~ObsError() {
  Log::trace() << "ObsError<OBS>::~ObsError starting" << std::endl;
  util::Timer timer(classname(), "~ObsError");
  err_.reset();
  Log::trace() << "ObsError<OBS>::~ObsError done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsError<OBS>::multiply(ObsVector_ & dy) const {
  Log::trace() << "ObsError<OBS>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  err_->multiply(dy);
  Log::trace() << "ObsError<OBS>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsError<OBS>::inverseMultiply(ObsVector_ & dy) const {
  Log::trace() << "ObsError<OBS>::inverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "inverseMultiply");
  err_->inverseMultiply(dy);
  Log::trace() << "ObsError<OBS>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsError<OBS>::sqrtMultiply(ObsVector_ & dy) const {
  Log::trace() << "ObsError<OBS>::sqrtMultiply starting" << std::endl;
  util::Timer timer(classname(), "sqrtMultiply");
  err_->sqrtMultiply(dy);
  Log::trace() << "ObsError<OBS>::sqrtMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsError<OBS>::invSqrtMultiply(ObsVector_ & dy) const {
  Log::trace() << "ObsError<OBS>::invSqrtMultiply starting" << std::endl;
  util::Timer timer(classname(), "invSqrtMultiply");
  err_->invSqrtMultiply(dy);
  Log::trace() << "ObsError<OBS>::invSqrtMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsError<OBS>::randomize(ObsVector_ & dy) const {
  Log::trace() << "ObsError<OBS>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  err_->randomize(dy);
  Log::trace() << "ObsError<OBS>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename OBS>
void ObsError<OBS>::save(const std::string & name) const {
  Log::trace() << "ObsError<OBS>::save starting" << std::endl;
  util::Timer timer(classname(), "save");
  err_->save(name);
  Log::trace() << "ObsError<OBS>::save done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
typename ObsError<OBS>::ObsVector_ ObsError<OBS>::obserrors() const {
  Log::trace() << "ObsError<OBS>::obserrors starting" << std::endl;
  util::Timer timer(classname(), "obserrors");
  ObsVector_ obserr = err_->obserrors();
  Log::trace() << "ObsError<OBS>::obserrors done" << std::endl;
  return obserr;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsError<OBS>::update(const ObsVector_ & obserr) {
  Log::trace() << "ObsError<OBS>::update starting" << std::endl;
  util::Timer timer(classname(), "update");
  err_->update(obserr);
  Log::trace() << "ObsError<OBS>::update done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsVector<OBS> ObsError<OBS>::inverseVariance() const {
  Log::trace() << "ObsError<OBS>::inverseVariance starting" << std::endl;
  util::Timer timer(classname(), "inverseVariance");
  ObsVector_ invar = err_->inverseVariance();
  Log::trace() << "ObsError<OBS>::inverseVariance done" << std::endl;
  return invar;
}

// -----------------------------------------------------------------------------

template <typename OBS>
double ObsError<OBS>::getRMSE() const {
  Log::trace() << "ObsError<OBS>::getRMSE starting" << std::endl;
  util::Timer timer(classname(), "getRMSE");
  double zz = err_->getRMSE();
  Log::trace() << "ObsError<OBS>::getRMSE done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsError<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsError<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << (*err_);
  Log::trace() << "ObsError<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERROR_H_
