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
#include "oops/interface/ObsAuxCovariance.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsAuxCovariances : public util::Printable,
                          private boost::noncopyable {
  typedef ObsAuxCovariance<MODEL>    ObsAuxCovariance_;
  typedef ObsAuxControls<MODEL>      ObsAuxControls_;
  typedef ObsAuxIncrements<MODEL>    ObsAuxIncrements_;

 public:
  static const std::string classname() {return "oops::ObsAuxCovariances";}

  explicit ObsAuxCovariances(const eckit::Configuration &);
  ~ObsAuxCovariances();

/// Operators
  void linearize(const ObsAuxControls_ &);
  void multiply(const ObsAuxIncrements_ &, ObsAuxIncrements_ &) const;
  void inverseMultiply(const ObsAuxIncrements_ &, ObsAuxIncrements_ &) const;
  void randomize(ObsAuxIncrements_ &) const;

  const eckit::LocalConfiguration & config() const {return conf_;}

 private:
  void print(std::ostream &) const;
  std::vector<boost::shared_ptr<ObsAuxCovariance_> > cov_;
  const eckit::LocalConfiguration conf_;
};

// =============================================================================

template<typename MODEL>
ObsAuxCovariances<MODEL>::ObsAuxCovariances(const eckit::Configuration & conf)
  : cov_(0), conf_(conf)
{
  Log::trace() << "ObsAuxCovariances<MODEL>::ObsAuxCovariances starting" << std::endl;
  std::vector<eckit::LocalConfiguration> obsconf;
  conf.get("ObsTypes", obsconf);
  for (std::size_t jobs = 0; jobs < obsconf.size(); ++jobs) {
    boost::shared_ptr<ObsAuxCovariance_>
        tmp(new ObsAuxCovariance_(obsconf[jobs]));
    cov_.push_back(tmp);
  }
  Log::trace() << "ObsAuxCovariances<MODEL>::ObsAuxCovariances done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsAuxCovariances<MODEL>::~ObsAuxCovariances() {
  Log::trace() << "ObsAuxCovariances<MODEL>::~ObsAuxCovariances starting" << std::endl;
  for (std::size_t jobs = 0; jobs < cov_.size(); ++jobs) {
    cov_[jobs].reset();
  }
  Log::trace() << "ObsAuxCovariances<MODEL>::~ObsAuxCovariances done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariances<MODEL>::linearize(const ObsAuxControls_ & xx) {
  Log::trace() << "ObsAuxCovariances<MODEL>::linearize starting" << std::endl;
  ASSERT(cov_.size() == xx.size());
  for (std::size_t jobs = 0; jobs < xx.size(); ++jobs) {
    cov_[jobs]->linearize(xx[jobs]);
  }
  Log::trace() << "ObsAuxCovariances<MODEL>::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariances<MODEL>::multiply(const ObsAuxIncrements_ & dx1,
                                       ObsAuxIncrements_ & dx2) const {
  Log::trace() << "ObsAuxCovariances<MODEL>::multiply starting" << std::endl;
  ASSERT(dx1.size() == dx2.size() && cov_.size() == dx1.size());
  for (std::size_t jobs = 0; jobs < dx1.size(); ++jobs) {
    cov_[jobs]->multiply(dx1[jobs], dx2[jobs]);
  }
  Log::trace() << "ObsAuxCovariances<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariances<MODEL>::inverseMultiply(const ObsAuxIncrements_ & dx1,
                                              ObsAuxIncrements_ & dx2) const {
  Log::trace() << "ObsAuxCovariances<MODEL>::inverseMultiply starting" << std::endl;
  ASSERT(dx1.size() == dx2.size() && cov_.size() == dx1.size());
  for (std::size_t jobs = 0; jobs < dx1.size(); ++jobs) {
    cov_[jobs]->inverseMultiply(dx1[jobs], dx2[jobs]);
  }
  Log::trace() << "ObsAuxCovariances<MODEL>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariances<MODEL>::randomize(ObsAuxIncrements_ & dx) const {
  Log::trace() << "ObsAuxCovariances<MODEL>::randomize starting" << std::endl;
  ASSERT(cov_.size() == dx.size());
  for (std::size_t jobs = 0; jobs < dx.size(); ++jobs) {
    cov_[jobs]->randomize(dx[jobs]);
  }
  Log::trace() << "ObsAuxCovariances<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariances<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxCovariances<MODEL>::print starting" << std::endl;
  for (std::size_t jobs = 0; jobs < cov_.size(); ++jobs) os << *cov_.at(jobs) << " ";
  Log::trace() << "ObsAuxCovariances<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSAUXCOVARIANCES_H_
