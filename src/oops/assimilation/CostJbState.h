/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTJBSTATE_H_
#define OOPS_ASSIMILATION_COSTJBSTATE_H_

#include <memory>
#include <vector>
#include <boost/noncopyable.hpp>

namespace oops {
  template<typename MODEL, typename OBS> class ControlVariable;
  template<typename MODEL, typename OBS> class ControlIncrement;
  template<typename MODEL> class Geometry;
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTLAD;
  template<typename MODEL> class PostProcessor;
  template<typename MODEL> class State;
  template<typename MODEL> class State4D;

// -----------------------------------------------------------------------------

/// Jb Cost Function Base Class
/*!
 * The CostJbState is the base class for the Jb term corresponding to the
 * state part of the control variable.
 */

template<typename MODEL, typename OBS> class CostJbState : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef Geometry<MODEL>               Geometry_;
  typedef State4D<MODEL>                State_;
  typedef JqTerm<MODEL>                 JqTerm_;
  typedef JqTermTLAD<MODEL>             JqTLAD_;

 public:
/// Constructor
  CostJbState() {}

/// Destructor
  virtual ~CostJbState() {}

  virtual void setPostProc(PostProcessor<State<MODEL>> &) {}
  virtual std::shared_ptr<JqTerm_> getJq() {return nullptr;}
  virtual void setTime(const CtrlVar_ &) {}

/// Get increment from state. This is usually first guess - background.
/// The third state argument is M(x) at the end of the window/subwindows for
/// computing the model error term (M(x_{i-1})-x_i) when active.
  virtual void computeIncrement(const CtrlVar_ &, const CtrlVar_ &, const std::shared_ptr<JqTerm_>,
                                CtrlInc_ &) const = 0;

/// Linearize before the linear computations.
  virtual void linearize(const CtrlVar_ &, const CtrlVar_ &, const Geometry_ &) = 0;

/// Add Jb gradient.
  virtual void addGradient(const CtrlInc_ &, CtrlInc_ &, CtrlInc_ &) const = 0;

/// Initialize Jq computations if needed.
  virtual JqTLAD_ * initializeJqTLAD() const = 0;

/// Finalize \f$ J_b\f$ after the TL run.
  virtual JqTLAD_ * initializeJqTL() const = 0;

/// Initialize \f$ J_b\f$ before the AD run.
  virtual JqTLAD_ * initializeJqAD(const CtrlInc_ &) const = 0;

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  virtual void Bmult(const CtrlInc_ &, CtrlInc_ &) const = 0;
  virtual void Bminv(const CtrlInc_ &, CtrlInc_ &) const = 0;

/// Randomize
  virtual void randomize(CtrlInc_ &) const = 0;

/// Accessors to data for constructing a new increment.
  virtual const Geometry_ & geometry() const = 0;
  virtual const Variables & variables() const = 0;
  virtual const std::vector<util::DateTime> & times() const = 0;
  virtual const eckit::mpi::Comm & comm() const = 0;
  virtual std::shared_ptr<State_> background() const = 0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBSTATE_H_
