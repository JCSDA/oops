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

#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

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
  typedef typename MODEL::Covariance Covariance_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const Variables &, const eckit::Configuration &,
                  const State_ &);
  virtual ~ErrorCovariance();

  void randomize(Increment_ &) const override;

 private:
  void doLinearize(const State_ &, const Geometry_ &) override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  void print(std::ostream &) const override;

  boost::scoped_ptr<Covariance_> covariance_;
};

// =============================================================================

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & resol, const Variables & vars,
                                        const eckit::Configuration & conf, const State_ & xb)
  : ModelSpaceCovarianceBase<MODEL>(conf), covariance_()
{
  Log::trace() << "ErrorCovariance<MODEL>::ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "ErrorCovariance");
  covariance_.reset(new Covariance_(resol.geometry(), vars, conf, xb.state()));
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
void ErrorCovariance<MODEL>::doLinearize(const State_ & xx, const Geometry_ & resol) {
  Log::trace() << "ErrorCovariance<MODEL>::doLinearize starting" << std::endl;
  util::Timer timer(classname(), "doLinearize");
  covariance_->linearize(xx.state(), resol.geometry());
  Log::trace() << "ErrorCovariance<MODEL>::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doMultiply(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "ErrorCovariance<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  covariance_->multiply(dx1.increment(), dx2.increment());
  Log::trace() << "ErrorCovariance<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doInverseMultiply(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");
  covariance_->inverseMultiply(dx1.increment(), dx2.increment());
  Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply done" << std::endl;
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
