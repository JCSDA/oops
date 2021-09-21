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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Parameters for diagonal obs errors
class ObsErrorDiagParameters : public ObsErrorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsErrorDiagParameters, ObsErrorParametersBase)
 public:
  /// perturbation amplitude multiplier
  Parameter<double> pert{"random amplitude", 1.0, this};
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
  ObsVector_ stddev_;
  ObsVector_ inverseVariance_;
  Parameters_ options_;
};

// =============================================================================

template<typename OBS>
ObsErrorDiag<OBS>::ObsErrorDiag(const ObsErrorDiagParameters & options, const ObsSpace_ & obsgeom)
  : stddev_(obsgeom, "ObsError"), inverseVariance_(obsgeom), options_(options)
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
