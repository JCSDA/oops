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

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/ObsErrorBase.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Diagonal observation error covariance matrix.

template<typename MODEL>
class ObsErrorDiag : public ObsErrorBase<MODEL> {
  typedef typename MODEL::ObsSpace              ObsSpace_;
  typedef typename MODEL::ObsVector             ObsVector_;

 public:
  ObsErrorDiag(const eckit::Configuration &, const ObsSpace_ &, const Variables &);
  ~ObsErrorDiag();

/// Multiply a Departure by \f$R\f$
  ObsVector_ * multiply(const ObsVector_ &) const;

/// Multiply a Departure by \f$R^{-1}\f$
  ObsVector_ * inverseMultiply(const ObsVector_ &) const;

/// Generate random perturbation
  void randomize(ObsVector_ &) const;

/// Get mean error for Jo table
  double getRMSE() const {return stddev_->rms();}

 private:
  void print(std::ostream &) const;

  boost::scoped_ptr<ObsVector_> stddev_;
  boost::scoped_ptr<ObsVector_> inverseVariance_;
};

// =============================================================================

template<typename MODEL>
ObsErrorDiag<MODEL>::ObsErrorDiag(const eckit::Configuration &, const ObsSpace_ & obsgeom,
                                  const Variables & observed)
  : ObsErrorBase<MODEL>(obsgeom, observed), stddev_(), inverseVariance_()
{
  stddev_.reset(new ObsVector_(obsgeom, observed));
  stddev_->read("EffectiveError");

  inverseVariance_.reset(new ObsVector_(*stddev_));
  *inverseVariance_ *= *inverseVariance_;
  inverseVariance_->invert();

  Log::debug() << "ObsErrorDiag:ObsErrorDiag constructed nobs = " << stddev_->nobs() << std::endl;
  Log::trace() << "ObsErrorDiag:ObsErrorDiag constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsErrorDiag<MODEL>::~ObsErrorDiag() {
  Log::trace() << "ObsErrorDiag:~ObsErrorDiag destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
typename MODEL::ObsVector * ObsErrorDiag<MODEL>::multiply(const ObsVector_ & dy) const {
  ObsVector_ * res = new ObsVector_(dy);
  *res /= *inverseVariance_;
  return res;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
typename MODEL::ObsVector * ObsErrorDiag<MODEL>::inverseMultiply(const ObsVector_ & dy) const {
  ObsVector_ * res = new ObsVector_(dy);
  *res *= *inverseVariance_;
  Log::debug() << "ObsErrorDiag:inverseMultiply nobs in = " << dy.nobs()
                               << ", out = " << res->nobs() << std::endl;
  return res;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorDiag<MODEL>::randomize(ObsVector_ & dy) const {
  dy.random();
  dy *= *stddev_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsErrorDiag<MODEL>::print(std::ostream & os) const {
  os << "ObsErrorDiag<MODEL>::print not implemeted yet";
}

// -----------------------------------------------------------------------------


}  // namespace oops

#endif  // OOPS_GENERIC_OBSERRORDIAG_H_
