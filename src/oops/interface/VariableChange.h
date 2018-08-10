/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_VARIABLECHANGE_H_
#define OOPS_INTERFACE_VARIABLECHANGE_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/base/VariableChangeBase.h"
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
class VariableChange : public oops::VariableChangeBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::VariableChange";}

  explicit VariableChange(const eckit::Configuration &);
  virtual ~VariableChange();

  void linearize(const State_ &, const Geometry_ &) override;

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
VariableChange<MODEL, CHVAR>::VariableChange(const eckit::Configuration & conf)
  : VariableChangeBase<MODEL>(conf), chvar_()
{
  Log::trace() << "VariableChange<MODEL, CHVAR>::VariableChange starting" << std::endl;
  util::Timer timer(classname(), "VariableChange");
  chvar_.reset(new CHVAR(conf));
  Log::trace() << "VariableChange<MODEL, CHVAR>::VariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
VariableChange<MODEL, CHVAR>::~VariableChange() {
  Log::trace() << "VariableChange<MODEL, CHVAR>::~VariableChange starting" << std::endl;
  util::Timer timer(classname(), "~VariableChange");
  chvar_.reset();
  Log::trace() << "VariableChange<MODEL, CHVAR>::~VariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void VariableChange<MODEL, CHVAR>::linearize(const State_ & xx, const Geometry_ & resol) {
  Log::trace() << "VariableChange<MODEL, CHVAR>::linearize starting" << std::endl;
  util::Timer timer(classname(), "linearize");
  chvar_->linearize(xx.state(), resol.geometry());
  Log::trace() << "VariableChange<MODEL, CHVAR>::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void VariableChange<MODEL, CHVAR>::multiply(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "VariableChange<MODEL, CHVAR>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  chvar_->multiply(dx1.increment(), dx2.increment());
  Log::trace() << "VariableChange<MODEL, CHVAR>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void VariableChange<MODEL, CHVAR>::multiplyInverse(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "VariableChange<MODEL, CHVAR>::multiplyInverse starting" << std::endl;
  util::Timer timer(classname(), "multiplyInverse");
  chvar_->multiplyInverse(dx1.increment(), dx2.increment());
  Log::trace() << "VariableChange<MODEL, CHVAR>::multiplyInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void VariableChange<MODEL, CHVAR>::multiplyAD(const Increment_ & dx1, Increment_ & dx2) const {
  Log::trace() << "VariableChange<MODEL, CHVAR>::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");
  chvar_->multiplyAD(dx1.increment(), dx2.increment());
  Log::trace() << "VariableChange<MODEL, CHVAR>::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void VariableChange<MODEL, CHVAR>::multiplyInverseAD(const Increment_ & dx1,
                                                     Increment_ & dx2) const {
  Log::trace() << "VariableChange<MODEL, CHVAR>::multiplyInverseAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyInverseAD");
  chvar_->multiplyInverseAD(dx1.increment(), dx2.increment());
  Log::trace() << "VariableChange<MODEL, CHVAR>::multiplyInverseAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void VariableChange<MODEL, CHVAR>::print(std::ostream & os) const {
  Log::trace() << "VariableChange<MODEL, CHVAR>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *chvar_;
  Log::trace() << "VariableChange<MODEL, CHVAR>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_VARIABLECHANGE_H_
