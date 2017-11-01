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

#include "util/Logger.h"
#include "oops/interface/ObsErrorBase.h"

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
  ObsErrorDiag(const ObsSpace_ &, const eckit::Configuration &);
  ~ObsErrorDiag();

/// Linearize and reset for inner loop (nothing in this simple case)
  void linearize(const ObsVector_ &) {}

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
ObsErrorDiag<MODEL>::ObsErrorDiag(const ObsSpace_ & obsgeom, const eckit::Configuration & config)
  : stddev_(), inverseVariance_()
{
  stddev_.reset(new ObsVector_(obsgeom));
  const std::string col = config.getString("obserror");
  stddev_->read(col);

  inverseVariance_.reset(new ObsVector_(*stddev_));
  *inverseVariance_ *= *inverseVariance_;
  inverseVariance_->invert();

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
