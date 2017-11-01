/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_ERRORCOVARIANCE_H_
#define OOPS_INTERFACE_ERRORCOVARIANCE_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

// Should factory be here and generic covariances wrtten at the MODEL::Increment level? YT

/// Wrapper for model space error covariances.
/*!
 *  This class provides the operations associated with the model space error
 *  covariance matrices (B or Q). It wraps the actual error covariance matrix
 *  which can be a model specific one or a generic one.
 */

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL> >,
                        private boost::noncopyable {
  typedef typename MODEL::Covariance        Covariance_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;
  typedef Variables<MODEL>           Variables_;

 public:
  static const std::string classname() {return "oops::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const Variables_ &, const eckit::Configuration &, const State_ &);
  virtual ~ErrorCovariance();

  void linearize(const State_ &, const Geometry_ &) override;
  void multiply(const Increment_ &, Increment_ &) const override;
  void inverseMultiply(const Increment_ &, Increment_ &) const override;
  void randomize(Increment_ &) const override;

 private:
  void print(std::ostream &) const override;
  boost::scoped_ptr<Covariance_> covariance_;
};

// =============================================================================

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & resol, const Variables_ & vars,
                                        const eckit::Configuration & conf, const State_ & xb)
  : covariance_()
{
  Log::trace() << "ErrorCovariance<MODEL>::ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "ErrorCovariance");
  covariance_.reset(new Covariance_(resol.geometry(), vars.variables(), conf, xb.state()));
  Log::trace() << "ErrorCovariance<MODEL>::ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::~ErrorCovariance() {
  Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance");
  covariance_.reset();
  Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::linearize(const State_ & xx, const Geometry_ & resol) {
  Log::trace() << "ErrorCovariance<MODEL>::linearize starting" << std::endl;
  util::Timer timer(classname(), "linearize");
  covariance_->linearize(xx.state(), resol.geometry());
  Log::trace() << "ErrorCovariance<MODEL>::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::multiply(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "ErrorCovariance<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  covariance_->multiply(dx1.increment(), dx2.increment());
  Log::trace() << "ErrorCovariance<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::inverseMultiply(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "ErrorCovariance<MODEL>::inverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "inverseMultiply");
  covariance_->inverseMultiply(dx1.increment(), dx2.increment());
  Log::trace() << "ErrorCovariance<MODEL>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::randomize(Increment_ & dx) const {
  Log::trace() << "ErrorCovariance<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");
  covariance_->randomize(dx.increment());
  Log::trace() << "ErrorCovariance<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ErrorCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *covariance_;
  Log::trace() << "ErrorCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_ERRORCOVARIANCE_H_
