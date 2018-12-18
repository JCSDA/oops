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
#include "oops/base/Variables.h"
#include "oops/interface/ObsErrorBase.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Diagonal observation error covariance matrix.

template<typename MODEL>
class ObsErrorDiag : public ObsErrorBase<MODEL> {
  typedef typename MODEL::ObsSpace              ObsSpace_;
  typedef typename MODEL::ObsVector             ObsVector_;

 public:
  ObsErrorDiag(const eckit::Configuration &, ObsSpace_ &, const Variables &);
  ~ObsErrorDiag();

/// Update after QC od other obs filters
  void update();

/// Multiply a Departure by \f$R\f$
  void multiply(ObsVector_ &) const;

/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(ObsVector_ &) const;

/// Generate random perturbation
  void randomize(ObsVector_ &) const;

/// Get mean error for Jo table
  double getRMSE() const {return stddev_.rms();}

 private:
  void print(std::ostream &) const;

  ObsVector_ stddev_;
  ObsVector_ inverseVariance_;
  double pert_;
};

// =============================================================================

template<typename MODEL>
ObsErrorDiag<MODEL>::ObsErrorDiag(const eckit::Configuration & conf, ObsSpace_ & obsgeom,
                                  const Variables & observed)
  : stddev_(obsgeom, observed), inverseVariance_(obsgeom, observed),
    pert_(conf.getDouble("random_amplitude", 1.0))
{
  stddev_.read("ObsError");
  stddev_.save("EffectiveError");

  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  Log::trace() << "ObsErrorDiag:ObsErrorDiag constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsErrorDiag<MODEL>::~ObsErrorDiag() {
  Log::trace() << "ObsErrorDiag:~ObsErrorDiag destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorDiag<MODEL>::update() {
  stddev_.read("EffectiveError");
  ObsVector_ qc(stddev_, false);
  qc.read("EffectiveQC");
  stddev_.mask(qc);
  stddev_.save("EffectiveError");

  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();

  Log::trace() << "ObsErrorDiag:update nobs = " << stddev_.nobs() << std::endl;
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
  dy *= pert_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorDiag<MODEL>::print(std::ostream & os) const {
  os << "ObsErrorDiag<MODEL>::print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_OBSERRORDIAG_H_
