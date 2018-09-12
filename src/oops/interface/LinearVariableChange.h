/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_LINEARVARIABLECHANGE_H_
#define OOPS_INTERFACE_LINEARVARIABLECHANGE_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/base/LinearVariableChangeBase.h"
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
/// Wrapper for change of variable

template <typename MODEL, typename CHVAR>
class LinearVariableChange : public oops::LinearVariableChangeBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::LinearVariableChange";}

  LinearVariableChange(const State_ &, const State_ &,
                       const Geometry_ &, const eckit::Configuration &);
  virtual ~LinearVariableChange();

  void multiply(const Increment_ &, Increment_ &) const override;
  void multiplyInverse(const Increment_ &, Increment_ &) const override;
  void multiplyAD(const Increment_ &, Increment_ &) const override;
  void multiplyInverseAD(const Increment_ &, Increment_ &) const override;

 private:
  void print(std::ostream &) const override;

  boost::scoped_ptr<CHVAR> chvar_;
};

// =============================================================================

template<typename MODEL, typename CHVAR>
LinearVariableChange<MODEL, CHVAR>::LinearVariableChange(const State_ & bg, const State_ & fg,
                                                         const Geometry_ & geom,
                                                         const eckit::Configuration & conf)
  : LinearVariableChangeBase<MODEL>(conf), chvar_()
{
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::LinearVariableChange starting" << std::endl;
  util::Timer timer(classname(), "LinearVariableChange");
  chvar_.reset(new CHVAR(bg.state(), fg.state(), geom.geometry(), conf));
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::LinearVariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
LinearVariableChange<MODEL, CHVAR>::~LinearVariableChange() {
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::~LinearVariableChange starting" << std::endl;
  util::Timer timer(classname(), "~LinearVariableChange");
  chvar_.reset();
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::~LinearVariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void LinearVariableChange<MODEL, CHVAR>::multiply(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  chvar_->multiply(dx1.increment(), dx2.increment());
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void LinearVariableChange<MODEL, CHVAR>::multiplyInverse(const Increment_ & dx1,
                                                         Increment_ & dx2) const {
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::multiplyInverse starting" << std::endl;
  util::Timer timer(classname(), "multiplyInverse");
  chvar_->multiplyInverse(dx1.increment(), dx2.increment());
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::multiplyInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void LinearVariableChange<MODEL, CHVAR>::multiplyAD(const Increment_ & dx1,
                                                    Increment_ & dx2) const {
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");
  chvar_->multiplyAD(dx1.increment(), dx2.increment());
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void LinearVariableChange<MODEL, CHVAR>::multiplyInverseAD(const Increment_ & dx1,
                                                           Increment_ & dx2) const {
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::multiplyInverseAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyInverseAD");
  chvar_->multiplyInverseAD(dx1.increment(), dx2.increment());
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::multiplyInverseAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void LinearVariableChange<MODEL, CHVAR>::print(std::ostream & os) const {
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *chvar_;
  Log::trace() << "LinearVariableChange<MODEL, CHVAR>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEARVARIABLECHANGE_H_
