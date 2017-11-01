/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSAUXCOVARIANCE_H_
#define OOPS_INTERFACE_OBSAUXCOVARIANCE_H_

#include <iostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "eckit/config/Configuration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsAuxCovariance : public util::Printable,
                         private boost::noncopyable,
                         private util::ObjectCounter<ObsAuxCovariance<MODEL> > {
  typedef typename MODEL::ObsAuxCovariance    ObsAuxCovariance_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncrement_;

 public:
  static const std::string classname() {return "oops::ObsAuxCovariance";}

  explicit ObsAuxCovariance(const eckit::Configuration &);
  ~ObsAuxCovariance();

/// Operators
  void linearize(const ObsAuxControl_ &);
  void multiply(const ObsAuxIncrement_ &, ObsAuxIncrement_ &) const;
  void inverseMultiply(const ObsAuxIncrement_ &, ObsAuxIncrement_ &) const;
  void randomize(ObsAuxIncrement_ &) const;

  const eckit::Configuration & config() const {return cov_->config();}

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<ObsAuxCovariance_> cov_;
};

// =============================================================================

template<typename MODEL>
ObsAuxCovariance<MODEL>::ObsAuxCovariance(const eckit::Configuration & conf) : cov_()
{
  Log::trace() << "ObsAuxCovariance<MODEL>::ObsAuxCovariance starting" << std::endl;
  util::Timer timer(classname(), "ObsAuxCovariance");
  cov_.reset(new ObsAuxCovariance_(conf));
  Log::trace() << "ObsAuxCovariance<MODEL>::ObsAuxCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsAuxCovariance<MODEL>::~ObsAuxCovariance() {
  Log::trace() << "ObsAuxCovariance<MODEL>::~ObsAuxCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ObsAuxCovariance");
  cov_.reset();
  Log::trace() << "ObsAuxCovariance<MODEL>::~ObsAuxCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariance<MODEL>::linearize(const ObsAuxControl_ & xx) {
  Log::trace() << "ObsAuxCovariance<MODEL>::linearize starting" << std::endl;
  util::Timer timer(classname(), "linearize");
  cov_->linearize(xx.obsauxcontrol());
  Log::trace() << "ObsAuxCovariance<MODEL>::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariance<MODEL>::multiply(const ObsAuxIncrement_ & dx1, ObsAuxIncrement_ & dx2) const {
  Log::trace() << "ObsAuxCovariance<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  cov_->multiply(dx1.obsauxincrement(), dx2.obsauxincrement());
  Log::trace() << "ObsAuxCovariance<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariance<MODEL>::inverseMultiply(const ObsAuxIncrement_ & dx1,
                                              ObsAuxIncrement_ & dx2) const {
  Log::trace() << "ObsAuxCovariance<MODEL>::inverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "inverseMultiply");
  cov_->inverseMultiply(dx1.obsauxincrement(), dx2.obsauxincrement());
  Log::trace() << "ObsAuxCovariance<MODEL>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariance<MODEL>::randomize(ObsAuxIncrement_ & dx) const {
  Log::trace() << "ObsAuxCovariance<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  cov_->randomize(dx.obsauxincrement());
  Log::trace() << "ObsAuxCovariance<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsAuxCovariance<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsAuxCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *cov_;
  Log::trace() << "ObsAuxCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSAUXCOVARIANCE_H_
