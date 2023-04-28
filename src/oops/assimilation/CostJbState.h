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
#include <boost/noncopyable.hpp>

namespace oops {

// Forward declaration
  template<typename MODEL> class Geometry;
  template<typename MODEL> class Increment;
  template<typename MODEL> class State;
  template<typename MODEL> class PostProcessor;
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

/// Jb Cost Function Base Class
/*!
 * The CostJbState is the base class for the Jb term corresponding to the
 * state part of the control variable.
 */

template<typename MODEL> class CostJbState : private boost::noncopyable {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;
  typedef PostProcessor<State_>      PostProc_;
  typedef JqTerm<MODEL>              JqTerm_;
  typedef JqTermTLAD<MODEL>          JqTLAD_;

 public:
/// Constructor
  CostJbState() {}

/// Destructor
  virtual ~CostJbState() {}

  virtual void setPostProc(PostProc_ &) {}
  virtual std::shared_ptr<JqTerm_> getJq() {return nullptr;}

/// Get increment from state. This is usually first guess - background.
/// The third state argument is M(x) at the end of the window/subwindows for
/// computing the model error term (M(x_{i-1})-x_i) when active.
  virtual void computeIncrement(const State_ &, const State_ &, const std::shared_ptr<JqTerm_>,
                                Increment_ &) const = 0;

/// Linearize before the linear computations.
  virtual void linearize(const State_ &, const Geometry_ &) = 0;

/// Add Jb gradient.
  virtual void addGradient(const Increment_ &, Increment_ &, Increment_ &) const = 0;

/// Initialize Jq computations if needed.
  virtual JqTLAD_ * initializeJqTLAD() const = 0;

/// Finalize \f$ J_b\f$ after the TL run.
  virtual JqTLAD_ * initializeJqTL() const = 0;

/// Initialize \f$ J_b\f$ before the AD run.
  virtual JqTLAD_ * initializeJqAD(const Increment_ &) const = 0;

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  virtual void Bmult(const Increment_ &, Increment_ &) const = 0;
  virtual void Bminv(const Increment_ &, Increment_ &) const = 0;

/// Randomize
  virtual void randomize(Increment_ &) const = 0;

/// Accessors to data for constructing a new increment.
  virtual const Geometry_ & geometry() const = 0;
  virtual const Variables & variables() const = 0;
  virtual const util::DateTime time() const = 0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBSTATE_H_
