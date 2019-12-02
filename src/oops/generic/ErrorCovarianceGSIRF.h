/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_ERRORCOVARIANCEGSIRF_H_
#define OOPS_GENERIC_ERRORCOVARIANCEGSIRF_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Timer.h"

namespace eckit {
  class LocalConfiguration;
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Model space error covariance on generic unstructured grid

template <typename MODEL>
class ErrorCovarianceGSIRF : public oops::ModelSpaceCovarianceBase<MODEL>,
                             private util::ObjectCounter<ErrorCovarianceGSIRF<MODEL>>,
                             private boost::noncopyable {
  typedef Geometry<MODEL>                         Geometry_;
  typedef Increment<MODEL>                        Increment_;
  typedef State<MODEL>                            State_;

 public:
  static const std::string classname() {return "oops::ErrorCovarianceGSIRF";}

  ErrorCovarianceGSIRF(const Geometry_ &, const Variables &,
                       const eckit::Configuration &, const State_ &, const State_ &);
  virtual ~ErrorCovarianceGSIRF();

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;
};

// =============================================================================

template<typename MODEL>
ErrorCovarianceGSIRF<MODEL>::ErrorCovarianceGSIRF(const Geometry_ & resol,
                                                  const Variables & vars,
                                                  const eckit::Configuration & conf,
                                                  const State_ & xb, const State_ & fg)
  : ModelSpaceCovarianceBase<MODEL>(xb, fg, resol, conf)
{
  Log::trace() << "ErrorCovarianceGSIRF::ErrorCovarianceGSIRF starting" << std::endl;

  Log::trace() << "ErrorCovarianceGSIRF::ErrorCovarianceGSIRF done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovarianceGSIRF<MODEL>::~ErrorCovarianceGSIRF() {
  Log::trace() << "ErrorCovarianceGSIRF<MODEL>::~ErrorCovarianceGSIRF starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovarianceGSIRF");
  Log::trace() << "ErrorCovarianceGSIRF<MODEL>::~ErrorCovarianceGSIRF done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceGSIRF<MODEL>::doRandomize(Increment_ & dx) const {
  Log::trace() << "ErrorCovarianceGSIRF<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");

  Log::trace() << "ErrorCovarianceGSIRF<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceGSIRF<MODEL>::doMultiply(const Increment_ & dxi,
                                             Increment_ & dxo) const {
  Log::trace() << "ErrorCovarianceGSIRF<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");

  dxo = dxi;

  Log::trace() << "ErrorCovarianceGSIRF<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceGSIRF<MODEL>::doInverseMultiply(const Increment_ & dxi,
                                                    Increment_ & dxo) const {
  Log::trace() << "ErrorCovarianceGSIRF<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");

  dxo = dxi;

  Log::trace() << "ErrorCovarianceGSIRF<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_ERRORCOVARIANCEGSIRF_H_
