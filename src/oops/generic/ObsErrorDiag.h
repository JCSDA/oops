/*
 * (C) Copyright 2009-2016 ECMWF.
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

#include "eckit/config/Configuration.h"
#include "oops/base/ObsErrorBase.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Parameters for diagonal obs errors
class ObsErrorDiagParameters : public Parameters {
 public:
  /// perturbation amplitude multiplier
  Parameter<double> pert{"random_amplitude", 1.0, this};
};

// -----------------------------------------------------------------------------
/// \brief Diagonal observation error covariance matrix.
template<typename MODEL>
class ObsErrorDiag : public ObsErrorBase<MODEL> {
  typedef ObsSpace<MODEL>              ObsSpace_;
  typedef ObsVector<MODEL>             ObsVector_;

 public:
  ObsErrorDiag(const eckit::Configuration &, const ObsSpace_ &);

/// Multiply a Departure by \f$R\f$
  void multiply(ObsVector_ &) const override;

/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(ObsVector_ &) const override;

/// Generate random perturbation
  void randomize(ObsVector_ &) const override;

/// Get mean error for Jo table
  double getRMSE() const override {return stddev_.rms();}

/// Return inverseVariance
  const ObsVector_ & inverseVariance() const {return inverseVariance_;}

 protected:
  ObsVector_ stddev_;
  ObsVector_ inverseVariance_;

 private:
  void print(std::ostream &) const override;
  ObsErrorDiagParameters options_;
};

// =============================================================================

template<typename MODEL>
ObsErrorDiag<MODEL>::ObsErrorDiag(const eckit::Configuration & conf, const ObsSpace_ & obsgeom)
  : stddev_(obsgeom, "EffectiveError"), inverseVariance_(obsgeom)
{
  options_.deserialize(conf);
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  Log::trace() << "ObsErrorDiag:ObsErrorDiag constructed nobs = " << stddev_.nobs() << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorDiag<MODEL>::multiply(ObsVector_ & dy) const {
  dy /= inverseVariance_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorDiag<MODEL>::inverseMultiply(ObsVector_ & dy) const {
  dy *= inverseVariance_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorDiag<MODEL>::randomize(ObsVector_ & dy) const {
  dy.random();
  dy *= stddev_;
  dy *= options_.pert;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorDiag<MODEL>::print(std::ostream & os) const {
  os << "Diagonal observation error covariance, inverse variances: "
     << inverseVariance_ << std::endl;
}

// -----------------------------------------------------------------------------


}  // namespace oops

#endif  // OOPS_GENERIC_OBSERRORDIAG_H_
