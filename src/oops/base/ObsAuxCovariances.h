/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSAUXCOVARIANCES_H_
#define OOPS_BASE_OBSAUXCOVARIANCES_H_

#include <iostream>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsAuxIncrements.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsAuxCovariance.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsAuxCovariances : public util::Printable,
                          private boost::noncopyable {
  typedef ObsAuxCovariance<OBS>    ObsAuxCovariance_;
  typedef ObsAuxControls<OBS>      ObsAuxControls_;
  typedef ObsAuxIncrements<OBS>    ObsAuxIncrements_;
  typedef ObsSpaces<OBS>           ObsSpaces_;

 public:
  static const std::string classname() {return "oops::ObsAuxCovariances";}

  ObsAuxCovariances(const ObsSpaces_ &, const eckit::Configuration &);
  ~ObsAuxCovariances();

/// Operators
  void linearize(const ObsAuxControls_ &);
  void multiply(const ObsAuxIncrements_ &, ObsAuxIncrements_ &) const;
  void inverseMultiply(const ObsAuxIncrements_ &, ObsAuxIncrements_ &) const;
  void randomize(ObsAuxIncrements_ &) const;

  const eckit::LocalConfiguration & config() const {return conf_;}
  const ObsSpaces_ & obspaces() const {return odb_;}

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsAuxCovariance_> > cov_;
  const ObsSpaces_ & odb_;
  const eckit::LocalConfiguration conf_;
};

// =============================================================================

template<typename OBS>
ObsAuxCovariances<OBS>::ObsAuxCovariances(const ObsSpaces_ & odb,
                                            const eckit::Configuration & conf)
  : cov_(0), odb_(odb), conf_(conf)
{
  Log::trace() << "ObsAuxCovariances<OBS>::ObsAuxCovariances starting" << std::endl;
  std::vector<eckit::LocalConfiguration> obsconf;
  conf.get("observations", obsconf);
  for (std::size_t jobs = 0; jobs < obsconf.size(); ++jobs) {
    boost::shared_ptr<ObsAuxCovariance_> tmp(new ObsAuxCovariance_(odb[jobs], obsconf[jobs]));
    cov_.push_back(tmp);
  }
  Log::trace() << "ObsAuxCovariances<OBS>::ObsAuxCovariances done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxCovariances<OBS>::~ObsAuxCovariances() {
  Log::trace() << "ObsAuxCovariances<OBS>::~ObsAuxCovariances starting" << std::endl;
  for (std::size_t jobs = 0; jobs < cov_.size(); ++jobs) {
    cov_[jobs].reset();
  }
  Log::trace() << "ObsAuxCovariances<OBS>::~ObsAuxCovariances done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxCovariances<OBS>::linearize(const ObsAuxControls_ & xx) {
  Log::trace() << "ObsAuxCovariances<OBS>::linearize starting" << std::endl;
  ASSERT(cov_.size() == xx.size());
  for (std::size_t jobs = 0; jobs < xx.size(); ++jobs) {
    cov_[jobs]->linearize(xx[jobs]);
  }
  Log::trace() << "ObsAuxCovariances<OBS>::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxCovariances<OBS>::multiply(const ObsAuxIncrements_ & dx1,
                                       ObsAuxIncrements_ & dx2) const {
  Log::trace() << "ObsAuxCovariances<OBS>::multiply starting" << std::endl;
  ASSERT(dx1.size() == dx2.size() && cov_.size() == dx1.size());
  for (std::size_t jobs = 0; jobs < dx1.size(); ++jobs) {
    cov_[jobs]->multiply(dx1[jobs], dx2[jobs]);
  }
  Log::trace() << "ObsAuxCovariances<OBS>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxCovariances<OBS>::inverseMultiply(const ObsAuxIncrements_ & dx1,
                                              ObsAuxIncrements_ & dx2) const {
  Log::trace() << "ObsAuxCovariances<OBS>::inverseMultiply starting" << std::endl;
  ASSERT(dx1.size() == dx2.size() && cov_.size() == dx1.size());
  for (std::size_t jobs = 0; jobs < dx1.size(); ++jobs) {
    cov_[jobs]->inverseMultiply(dx1[jobs], dx2[jobs]);
  }
  Log::trace() << "ObsAuxCovariances<OBS>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxCovariances<OBS>::randomize(ObsAuxIncrements_ & dx) const {
  Log::trace() << "ObsAuxCovariances<OBS>::randomize starting" << std::endl;
  ASSERT(cov_.size() == dx.size());
  for (std::size_t jobs = 0; jobs < dx.size(); ++jobs) {
    cov_[jobs]->randomize(dx[jobs]);
  }
  Log::trace() << "ObsAuxCovariances<OBS>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxCovariances<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxCovariances<OBS>::print starting" << std::endl;
  for (std::size_t jobs = 0; jobs < cov_.size(); ++jobs) os << *cov_.at(jobs) << " ";
  Log::trace() << "ObsAuxCovariances<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSAUXCOVARIANCES_H_
