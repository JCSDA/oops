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

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Departures.h"
#include "oops/base/ObsError.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Container for ObsErrors for all observation types that are used in DA
template <typename OBS>
class ObsErrors : public util::Printable,
                  private boost::noncopyable {
  typedef Departures<OBS>          Departures_;
  typedef ObsError<OBS>            ObsError_;
  typedef ObsSpaces<OBS>           ObsSpaces_;

 public:
  static const std::string classname() {return "oops::ObsErrors";}

  ObsErrors(const eckit::Configuration &, const ObsSpaces_ &);

/// Accessor and size
  size_t size() const {return err_.size();}
  ObsError_ & operator[](const size_t ii) {return err_.at(ii);}
  const ObsError_ & operator[](const size_t ii) const {return err_.at(ii);}

/// Multiply a Departure by \f$R\f$
  void multiply(Departures_ &) const;
/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(Departures_ &) const;

/// Generate random perturbation
  void randomize(Departures_ &) const;

/// Save obs errors
  void save(const std::string &) const;

  /// returns inverse of observation error variance
  Departures_ inverseVariance() const;

 private:
  void print(std::ostream &) const override;
  std::vector<ObsError_> err_;
  const ObsSpaces_ & os_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsErrors<OBS>::ObsErrors(const eckit::Configuration & config,
                          const ObsSpaces_ & os) : err_(), os_(os) {
  std::vector<eckit::LocalConfiguration> obsconf = config.getSubConfigurations();
  for (size_t jj = 0; jj < os.size(); ++jj) {
    eckit::LocalConfiguration conf = obsconf[jj].getSubConfiguration("obs error");
    err_.emplace_back(conf, os[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsErrors<OBS>::multiply(Departures_ & dy) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj].multiply(dy[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsErrors<OBS>::inverseMultiply(Departures_ & dy) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj].inverseMultiply(dy[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsErrors<OBS>::randomize(Departures_ & dy) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj].randomize(dy[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsErrors<OBS>::save(const std::string & name) const {
  for (const auto & err : err_) {
    err.save(name);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
Departures<OBS> ObsErrors<OBS>::inverseVariance() const {
  Departures_ invvar(os_);
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    invvar[jj] = err_[jj].inverseVariance();
  }
  return invvar;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrors<OBS>::print(std::ostream & os) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) os << err_[jj] << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERRORS_H_
