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

#include "oops/base/Geometry.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"

namespace eckit {
  class Configuration;
}

namespace oops {

template <typename MODEL, typename OBS> class ControlIncrement;
template <typename MODEL, typename OBS> class ControlVariable;

// -----------------------------------------------------------------------------

/// Base Class for Cost Function Terms
/*!
 * Abstract base class for the terms of the cost function (other than Jb).
 */

template<typename MODEL, typename OBS> class CostTermBase {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef PostProcessor<State_>      PostProc_;
  typedef PostProcessorTLAD<MODEL>   PostProcTLAD_;

 public:
/// Destructor
  virtual ~CostTermBase() {}

/// Initialize and set post-processors to collect data during nonlinear model integration
  virtual void setPostProc(const ControlVariable<MODEL, OBS> &, const eckit::Configuration &,
                           PostProc_ &) = 0;
/// Finish computation of cost function term after nonlinear model integration
  virtual double computeCost() = 0;
  virtual void printCostTestHack() = 0;

/// Set post-processors for nonlinear model integration and save linearisation trajectory
  virtual void setPostProcTraj(const ControlVariable<MODEL, OBS> &, const eckit::Configuration &,
                               const Geometry_ &, PostProcTLAD_ &) = 0;
/// Finish cost computation and trajectory handling after nonlinear model integration
  virtual void computeCostTraj() = 0;

/// Initialize and set TL post-processors to collect data during TL model integration
  virtual void setPostProcTL(const ControlIncrement<MODEL, OBS> &, PostProcTLAD_ &) const = 0;
/// Finish cost computation after TL model integration
  virtual void computeCostTL(const ControlIncrement<MODEL, OBS> &,
                             GeneralizedDepartures &) const = 0;

/// Adjoint of computeCostTL (initialize and set post-processors adjoint to force AD model)
// Going by the book, computeCostAD should have the same arguments as computeCostTL (with
// swapped constness). PostProcTLAD is added because it is the way forcing is passed to the
// model (adjoint operations are called in reverse order so computeCostAD will come first).
  virtual void computeCostAD(std::shared_ptr<const GeneralizedDepartures>,
                             ControlIncrement<MODEL, OBS> &, PostProcTLAD_ &) const = 0;
/// Adjoint ot setPostProcTL (clean-up)
  virtual void setPostProcAD() const = 0;

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
