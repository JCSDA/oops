/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTTERMBASE_H_
#define OOPS_ASSIMILATION_COSTTERMBASE_H_

#include <memory>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Base Class for Cost Function Terms
/*!
 * Abstract base class for the terms of the cost function.
 */

template<typename MODEL, typename OBS> class CostTermBase {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Increment<MODEL>           Increment_;
  typedef PostProcessor<State_>      PostProc_;
  typedef PostProcessorTLAD<MODEL>   PostProcTLAD_;

 public:
/// Destructor
  virtual ~CostTermBase() {}

/// Initialize before nonlinear model integration.
  virtual void initialize(const ControlVariable<MODEL, OBS> &, const eckit::Configuration &,
                          PostProc_ &) = 0;
  virtual void initializeTraj(const ControlVariable<MODEL, OBS> &, const Geometry_ &,
                              const eckit::Configuration &, PostProcTLAD_ &) = 0;

/// Finalize computation after nonlinear model integration.
  virtual double finalize() = 0;
  virtual void finalizeTraj() = 0;

/// Initialize before starting the TL run.
  virtual void setupTL(const ControlIncrement<MODEL, OBS> &, PostProcTLAD_ &) const = 0;

/// Initialize before starting the AD run.
  virtual void setupAD(std::shared_ptr<const GeneralizedDepartures>,
                       ControlIncrement<MODEL, OBS> &, PostProcTLAD_ &) const = 0;

/// Multiply by covariance (or weight) matrix and its inverse.
  virtual std::unique_ptr<GeneralizedDepartures>
    multiplyCovar(const GeneralizedDepartures &) const = 0;
  virtual std::unique_ptr<GeneralizedDepartures>
    multiplyCoInv(const GeneralizedDepartures &) const = 0;

/// Provide new dual space vector (for example a Departure for Jo).
  virtual std::unique_ptr<GeneralizedDepartures> newDualVector() const = 0;

/// Gradient at first guess.
  virtual std::unique_ptr<GeneralizedDepartures> newGradientFG() const = 0;

/// Reset trajectory.
  virtual void resetLinearization() = 0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTTERMBASE_H_
