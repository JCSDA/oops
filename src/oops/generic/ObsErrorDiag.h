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
#include <string>

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
  OOPS_CONCRETE_PARAMETERS(ObsErrorDiagParameters, Parameters)
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
  ObsErrorDiag(const eckit::Configuration &, const ObsSpace_ &);

/// Update after obs errors potentially changed
  void update() override;

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
  ObsVector_ & obserrors() override {return stddev_;}
  const ObsVector_ & obserrors() const override {return stddev_;}

/// Return inverseVariance
  const ObsVector_ & inverseVariance() const override {return inverseVariance_;}

 protected:
  ObsVector_ stddev_;
  ObsVector_ inverseVariance_;

 private:
  void print(std::ostream &) const override;
  ObsErrorDiagParameters options_;
};

// =============================================================================

template<typename OBS>
ObsErrorDiag<OBS>::ObsErrorDiag(const eckit::Configuration & conf, const ObsSpace_ & obsgeom)
  : stddev_(obsgeom, "ObsError"), inverseVariance_(obsgeom)
{
  options_.deserialize(conf);
  this->update();
  Log::trace() << "ObsErrorDiag:ObsErrorDiag constructed nobs = " << stddev_.nobs() << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrorDiag<OBS>::update() {
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
