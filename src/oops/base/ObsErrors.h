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
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsError.h"
#include "oops/util/ConfigFunctions.h"  // for vectoriseAndFilter
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Container for ObsErrors for all observation types that are used in DA
template <typename OBS>
class ObsErrors : public util::Printable,
                  private boost::noncopyable {
  typedef Departures<OBS>                Departures_;
  typedef ObsError<OBS>                  ObsError_;
  typedef ObsSpaces<OBS>                 ObsSpaces_;

 public:
  static const std::string classname() {return "oops::ObsErrors";}

  ObsErrors(const std::vector<eckit::LocalConfiguration> &, const ObsSpaces_ &);
  ObsErrors(const eckit::Configuration &, const ObsSpaces_ &);

/// Accessor and size
  size_t size() const {return err_.size();}
  ObsError_ & operator[](const size_t ii) {return err_.at(ii);}
  const ObsError_ & operator[](const size_t ii) const {return err_.at(ii);}

/// Multiply a Departure by \f$R\f$
  void multiply(Departures_ &) const;
/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(Departures_ &) const;

  Eigen::MatrixXd localInverseMultiply(const Eigen::MatrixXf &) const;

/// Create local R matrix
  void localize(Departures_ &) const;

  /// Generate random perturbation
  void randomize(Departures_ &) const;

/// Save obs errors
  void save(const std::string &) const;

  /// returns inverse of observation error variance
  Departures_ inverseVariance() const;

  Eigen::VectorXd local_invVarR() const;

 private:
  void print(std::ostream &) const override;
  std::vector<ObsError_> err_;
  const ObsSpaces_ & os_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsErrors<OBS>::ObsErrors(const std::vector<eckit::LocalConfiguration> & obsConfs,
                          const ObsSpaces_ & os) : err_(), os_(os) {
  ASSERT(obsConfs.empty() || obsConfs.size() == os.size());
  const eckit::LocalConfiguration defaultObsConf;

  err_.reserve(os.size());
  for (size_t jj = 0; jj < os.size(); ++jj) {
    const eckit::LocalConfiguration & obsConf = obsConfs.empty() ? defaultObsConf : obsConfs[jj];
    err_.emplace_back(obsConf, os_[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsErrors<OBS>::ObsErrors(const eckit::Configuration & config, const ObsSpaces_ & os)
  : err_(), os_(os)
{
  std::vector<eckit::LocalConfiguration> subconfigs = config.getSubConfigurations();
  ASSERT(subconfigs.size() == os.size());
  err_.reserve(os.size());
  for (size_t jj = 0; jj < os.size(); ++jj) {
    eckit::LocalConfiguration obsConf;
    if (subconfigs[jj].has("obs error")) {
      obsConf = subconfigs[jj].getSubConfiguration("obs error");
    }
    err_.emplace_back(obsConf, os_[jj]);
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
Eigen::MatrixXd ObsErrors<OBS>::localInverseMultiply(const Eigen::MatrixXf & zz) const {
  size_t i_obs = 0;
  Eigen::MatrixXf zz_jj;
  Eigen::MatrixXf zzRinv(zz.rows(), zz.cols());
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    // get components of zz from this obs space.
    // Assumption is that ordering of obs spaces is consistent between
    // zz and this object.
    zz_jj = zz.block(0, i_obs, zz.rows(), err_[jj].localDim());
    zzRinv.block(0, i_obs, zz.rows(), err_[jj].localDim()) =
        err_[jj].localInverseMultiply(zz_jj);
    i_obs += err_[jj].localDim();
  }

  return zzRinv.cast<double>();
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsErrors<OBS>::localize(Departures_ & locvector) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    err_[jj].localize(locvector[jj]);
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

template <typename OBS>
Eigen::VectorXd ObsErrors<OBS>::local_invVarR() const {
  size_t size = 0;
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    size += err_[jj].local_invVarR().size();
  }
  Eigen::VectorXd local_invVar(size);
  size_t offset = 0;
  for (size_t jj = 0; jj < err_.size(); ++jj) {
    size_t size_jj = err_[jj].local_invVarR().size();
    local_invVar.middleRows(offset, size_jj) = err_[jj].local_invVarR();
    offset += size_jj;
  }
  return local_invVar;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsErrors<OBS>::print(std::ostream & os) const {
  for (size_t jj = 0; jj < err_.size(); ++jj) os << err_[jj] << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERRORS_H_
