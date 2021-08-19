/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSERRORBASE_H_
#define OOPS_INTERFACE_OBSERRORBASE_H_

#include <memory>
#include <string>

#include "eckit/system/ResourceUsage.h"

#include "oops/generic/ObsErrorBase.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

namespace interface {

// -----------------------------------------------------------------------------
/// \brief Base class for OBS-specific implementations of the ObsError interface.
///
/// interface::ObsErrorBase overrides oops::ObsErrorBase methods to pass OBS-specific
/// implementations of ObsVector to OBS-specific implementations of ObsError.
///
/// Note: each subclass should typedef `Parameters_` to the name of a subclass of
/// ObsErrorParametersBase holding its configuration settings and provide a constructor with the
/// following signature:
///
///     ObsErrorBase(const Parameters_ &params, OBS::ObsSpace &obsspace,
///                  const eckit::mpi::Comm &timeComm);
///
/// This constructor should pass \c timeComm to the constructor of this class.
template <typename OBS>
class ObsErrorBase : public oops::ObsErrorBase<OBS> {
  typedef typename OBS::ObsSpace  ObsSpace_;
  typedef typename OBS::ObsVector ObsVector_;

 public:
  explicit ObsErrorBase(const eckit::mpi::Comm &timeComm)
    : timeComm_(timeComm) {}

  // Overrides of oops::ObsErrorBase methods, converting between oops interface classes and
  // OBS-specific classes operated upon by OBS-specific implementations of the ObsError interface.

  void multiply(ObsVector<OBS> &dy) const final {
    this->multiply(dy.obsvector());
  }

  void inverseMultiply(ObsVector<OBS> &dy) const final {
    this->inverseMultiply(dy.obsvector());
  }

  void randomize(ObsVector<OBS> &dy) const final {
    this->randomize(dy.obsvector());
  }

  ObsVector<OBS> obserrors() const final {
    return ObsVector<OBS>(this->getObsErrors(), timeComm_);
  }

  void update(const ObsVector<OBS> &dy) final {
    this->update(dy.obsvector());
  }

  ObsVector<OBS> inverseVariance() const final {
    return ObsVector<OBS>(this->getInverseVariance(), timeComm_);
  }

  // The methods below need to be overridden in subclasses (along with save(), getRMSE() and
  // print(), methods inherited from parent classes that neither take nor return instances of oops
  // interface classes and therefore are still abstract).

  /// Multiply a Departure \p dy by \f$R\f$.
  virtual void multiply(ObsVector_ &dy) const = 0;

  /// Multiply a Departure \p dy by \f$R^{-1}\f$.
  virtual void inverseMultiply(ObsVector_ &dy) const = 0;

  /// Generate a random perturbation in \p dy.
  virtual void randomize(ObsVector_ &dy) const = 0;

  /// Return a copy of obs error std. dev. If this ObsVector_ is modified (e.g. by obs filters),
  /// it should be passed back to update() to ensure the covariance matrix stays consistent.
  virtual std::unique_ptr<ObsVector_> getObsErrors() const = 0;

  /// Set the diagonal of the covariance matrix to \p stddev squared.
  virtual void update(const ObsVector_ &stddev) = 0;

  /// Return the vector of inverse obs error variances.
  virtual std::unique_ptr<ObsVector_> getInverseVariance() const = 0;

 private:
  const eckit::mpi::Comm & timeComm_;
};

// =============================================================================

/// \brief A subclass of ObsErrorFactory able to create instances of T (a concrete subclass of
/// interface::ObsErrorBase<OBS>). Passes OBS::ObsSpace to the constructor of T.
template <class OBS, class T>
class ObsErrorMaker : public ObsErrorFactory<OBS> {
 public:
  typedef oops::ObsErrorBase<OBS>    ObsErrorBase_;
  typedef oops::ObsErrorFactory<OBS> ObsErrorFactory_;
  typedef oops::ObsSpace<OBS>        ObsSpace_;
  typedef typename T::Parameters_    Parameters_;

  explicit ObsErrorMaker(const std::string & name) : ObsErrorFactory_(name) {}

 private:
  std::unique_ptr<ObsErrorBase_> make(const ObsErrorParametersBase & parameters,
                                      const ObsSpace_ & obsspace) override {
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return std::make_unique<T>(stronglyTypedParameters, obsspace.obsspace(), obsspace.timeComm());
  }

  std::unique_ptr<ObsErrorParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }
};

// -----------------------------------------------------------------------------

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSERRORBASE_H_
