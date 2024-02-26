/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_OBSERRORDIAG_H_
#define OOPS_GENERIC_OBSERRORDIAG_H_

#include <sstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/base/ObsVector.h"
#include "oops/generic/ObsErrorBase.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Parameters for diagonal obs errors
class ObsErrorDiagParameters : public ObsErrorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsErrorDiagParameters, ObsErrorParametersBase)

 public:
  /// perturbation amplitude multiplier
  Parameter<double> pert{"obs perturbations amplitude", 1.0, this};

  /// Set to true to constrain observation perturbations to have a zero ensemble mean.
  ///
  /// Important: for this to work, the following requirements must be satisfied:
  ///
  /// 1. The `obs perturbations seed` option in the `obs space` section must be set to the same
  /// value for all ensemble members.
  ///
  /// 2. All ensemble members must use the same observations in the same order.
  Parameter<bool> zeroMeanPerturbations{"zero-mean perturbations", false, this};

  /// 1-based ensemble member index.
  ///
  /// Used (and required) only if `zero-mean perturbations` is set to true.
  OptionalParameter<int> member{"member", this, {minConstraint(1)}};

  /// Number of ensemble members.
  ///
  /// Used (and required) only if `zero-mean perturbations` is set to true.
  OptionalParameter<int> numberOfMembers{"number of members", this, {minConstraint(1)}};

  // Import both overloads from the base class; we're going to override one of them
  using ObsErrorParametersBase::deserialize;

  /// Overridden to detect missing conditionally required parameters
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;
};

// -----------------------------------------------------------------------------
/// \brief Diagonal observation error covariance matrix.
template<typename OBS>
class ObsErrorDiag : public ObsErrorBase<OBS> {
  typedef ObsSpace<OBS>              ObsSpace_;
  typedef ObsVector<OBS>             ObsVector_;

 public:
  /// The type of parameters passed to the constructor.
  /// This typedef is used by the ObsErrorFactory.
  typedef ObsErrorDiagParameters Parameters_;

  ObsErrorDiag(const Parameters_ &, const ObsSpace_ &);

/// Update after obs errors potentially changed
  void update(const ObsVector_ &) override;

/// Multiply a Departure by \f$R\f$
  void multiply(ObsVector_ &) const override;

/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(ObsVector_ &) const override;

/// Generate random perturbation
  void randomize(ObsVector_ &) const override;

/// Save obs errors
  void save(const std::string &) const override;

/// Get mean std deviation of errors for Jo table
  double getRMSE() const override {return stddev_.rms();}

/// Get obs errors std deviation
  ObsVector_ obserrors() const override {return stddev_;}

/// Get inverseVariance
  ObsVector_ inverseVariance() const override {return inverseVariance_;}

 private:
  void print(std::ostream &) const override;

  void randomizeWithoutZeroEnsembleMean(ObsVector_ &) const;

  void randomizeWithZeroEnsembleMean(ObsVector_ &) const;

  ObsVector_ stddev_;
  ObsVector_ inverseVariance_;
  Parameters_ options_;
};

// =============================================================================

template<typename OBS>
ObsErrorDiag<OBS>::ObsErrorDiag(const ObsErrorDiagParameters & options, const ObsSpace_ & obsgeom)
  : stddev_(obsgeom, "ObsError"), inverseVariance_(obsgeom, ""), options_(options)
{
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  Log::trace() << "ObsErrorDiag:ObsErrorDiag constructed nobs = " << stddev_.nobs() << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::update(const ObsVector_ & obserr) {
  stddev_ = obserr;
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  Log::info() << "ObsErrorDiag covariance updated " << stddev_.nobs() << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::multiply(ObsVector_ & dy) const {
  dy /= inverseVariance_;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::inverseMultiply(ObsVector_ & dy) const {
  dy *= inverseVariance_;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::randomize(ObsVector_ & dy) const {
  if (options_.zeroMeanPerturbations)
    randomizeWithZeroEnsembleMean(dy);
  else
    randomizeWithoutZeroEnsembleMean(dy);
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::randomizeWithZeroEnsembleMean(ObsVector_ & dy) const {
  ObsVector_ perturbation(dy);
  ObsVector_ sum(dy);
  sum.zero();

  // Generate initial independent perturbations for all ensemble members.
  // Calculate their sum and store this member's perturbations in 'dy'.
  const int myMember = options_.member.value().value();
  const int numMembers = options_.numberOfMembers.value().value();
  for (int member = 1; member <= numMembers; ++member) {
    perturbation.random();
    sum += perturbation;
    if (member == myMember)
      dy = perturbation;
  }

  // Subtract the ensemble mean of perturbations from this member's perturbations.
  dy.axpy(-1.0 / numMembers, sum);

  // Scale perturbations to the requested amplitude.
  dy *= stddev_;
  dy *= std::sqrt(numMembers / (numMembers - 1.0)) * options_.pert.value();
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::randomizeWithoutZeroEnsembleMean(ObsVector_ & dy) const {
  dy.random();
  dy *= stddev_;
  dy *= options_.pert;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::save(const std::string & name) const {
  stddev_.save(name);
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::print(std::ostream & os) const {
  os << "Diagonal observation error covariance" << std::endl << stddev_;
}

// -----------------------------------------------------------------------------


}  // namespace oops

#endif  // OOPS_GENERIC_OBSERRORDIAG_H_
