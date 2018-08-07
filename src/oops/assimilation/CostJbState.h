/*
 * (C) Copyright 2009-2016 ECMWF.
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

#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"

namespace util {
  class Duration;
}

namespace oops {

// Forward declaration
  template<typename MODEL> class ControlIncrement;
  template<typename MODEL> class Increment4D;
  template<typename MODEL> class State4D;
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

/// Jb Cost Function Base Class
/*!
 * The CostJbState is the base class for the Jb term corresponding to the
 * state part (3D or 4D) of the control variable.
 */

template<typename MODEL> class CostJbState : private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef State4D<MODEL>             State4D_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef Geometry<MODEL>            Geometry_;

 public:
/// Constructor
  CostJbState() {}

/// Destructor
  virtual ~CostJbState() {}

/// Initialize Jq computations if needed.
  virtual JqTerm<MODEL> * initializeJq() const = 0;
  virtual JqTermTLAD<MODEL> * initializeJqTLAD() const = 0;

/// Get increment from state (usually first guess).
  virtual void computeIncrement(const State4D_ &, const State4D_ &,
                                Increment4D_ &) const = 0;

/// Linearize before the linear computations.
  virtual void linearize(const State4D_ &, const Geometry_ &) = 0;

/// Add Jb gradient.
  virtual void addGradient(const Increment4D_ &, Increment4D_ &, Increment4D_ &) const = 0;

/// Finalize \f$ J_b\f$ after the TL run.
  virtual JqTermTLAD<MODEL> * initializeJqTL() const = 0;

/// Initialize \f$ J_b\f$ before the AD run.
  virtual JqTermTLAD<MODEL> * initializeJqAD(const Increment4D_ &) const = 0;

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  virtual void Bmult(const Increment4D_ &, Increment4D_ &) const = 0;
  virtual void Bminv(const Increment4D_ &, Increment4D_ &) const = 0;

/// Multiply by \f$ K\f$ and \f$ K^T\f$  
  virtual Increment4D_ Kmult(const Increment4D_ &) const = 0;
  virtual Increment4D_ KmultAdjoint(const Increment4D_ &) const = 0;

/// Randomize
  virtual void randomize(Increment4D_ &) const = 0;

/// Create new increment (set to 0).
  virtual unsigned int nstates() const = 0;
  virtual Increment_ * newStateIncrement(const unsigned int) const = 0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBSTATE_H_
