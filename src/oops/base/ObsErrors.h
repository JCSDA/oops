/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERRORS_H_
#define OOPS_BASE_OBSERRORS_H_

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Departures.h"
#include "oops/base/ObsErrorBase.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \biref Container for ObsErrors for all observation types that are used in DA
template <typename OBS>
class ObsErrors : public util::Printable,
                  private boost::noncopyable {
  typedef Departures<OBS>          Departures_;
  typedef ObsErrorBase<OBS>        ObsError_;
  typedef ObsSpaces<OBS>           ObsSpaces_;
  typedef ObsVector<OBS>           ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsErrors";}

  ObsErrors(const eckit::Configuration &, const ObsSpaces_ &);

/// Accessor and size
  size_t size() const {return err_.size();}
  const ObsError_ & operator[](const size_t ii) const {return *err_.at(ii);}

/// Multiply a Departure by \f$R\f$
  void multiply(Departures_ &) const;
/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(Departures_ &) const;

/// Generate random perturbation
  void randomize(Departures_ &) const;

/// Pack inverseVariance into an Eigen vector (excluding observations
///  that are masked out)
  Eigen::VectorXd packInverseVarianceEigen() const;

 private:
  void print(std::ostream &) const;
  std::vector<std::unique_ptr<ObsError_> > err_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsErrors<OBS>::ObsErrors(const eckit::Configuration & config,
                            const ObsSpaces_ & os) : err_() {
  std::vector<eckit::LocalConfiguration> obsconf;
  config.get("observations", obsconf);
  for (size_t jj = 0; jj < os.size(); ++jj) {
    eckit::LocalConfiguration conf(obsconf[jj], "obs error");
    err_.emplace_back(ObsErrorFactory<OBS>::create(conf, os[jj]));
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsErrors<OBS>::multiply(Departures_ & dy) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj]->multiply(dy[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsErrors<OBS>::inverseMultiply(Departures_ & dy) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj]->inverseMultiply(dy[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsErrors<OBS>::randomize(Departures_ & dy) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj]->randomize(dy[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
Eigen::VectorXd ObsErrors<OBS>::packInverseVarianceEigen() const {
  // compute nobs accross all obs errors
  unsigned int nobs = 0;
  for (size_t iov = 0; iov < err_.size(); ++iov) {
    const ObsVector_ & ov = err_[iov]->inverseVariance();
    nobs += ov.nobs();
  }

  // concatinate all inverseVariance into a 1d vector
  Eigen::VectorXd vec(nobs);
  unsigned int ii = 0;
  for (size_t iov = 0; iov < err_.size(); ++iov) {
    const ObsVector_ & ov = err_[iov]->inverseVariance();
    vec.segment(ii, ov.nobs()) = ov.packEigen();
    ii += ov.nobs();
  }
  ASSERT(ii == nobs);
  return vec;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrors<OBS>::print(std::ostream & os) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) os << *err_[jj];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERRORS_H_
