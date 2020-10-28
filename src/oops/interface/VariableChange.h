/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_VARIABLECHANGE_H_
#define OOPS_INTERFACE_VARIABLECHANGE_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/VariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
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
  typedef State<MODEL>               State_;

 public:
  /// Defined as CHVAR::Parameters_ if CHVAR defines a Parameters_ type; otherwise as
  /// GenericVariableChangeParameters
  typedef TParameters_IfAvailableElseFallbackType_t<
    CHVAR, GenericVariableChangeParameters> Parameters_;

  static const std::string classname() {return "oops::VariableChange";}

  VariableChange(const Geometry_ &, const Parameters_ &);
  virtual ~VariableChange();

  void changeVar(const State_ &, State_ &) const override;
  void changeVarInverse(const State_ &, State_ &) const override;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<CHVAR> chvar_;
};

// =============================================================================

template<typename MODEL, typename CHVAR>
VariableChange<MODEL, CHVAR>::VariableChange(const Geometry_ & geom,
                                             const Parameters_ & params)
  : VariableChangeBase<MODEL>(params), chvar_()
{
  Log::trace() << "VariableChange<MODEL, CHVAR>::VariableChange starting" << std::endl;
  util::Timer timer(classname(), "VariableChange");
  chvar_.reset(new CHVAR(geom.geometry(),
                         parametersOrConfiguration<HasParameters_<CHVAR>::value>(params)));
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
void VariableChange<MODEL, CHVAR>::changeVar(const State_ & x1, State_ & x2) const {
  Log::trace() << "VariableChange<MODEL, CHVAR>::changeVar starting" << std::endl;
  util::Timer timer(classname(), "changeVar");
  chvar_->changeVar(x1.state(), x2.state());
  Log::trace() << "VariableChange<MODEL, CHVAR>::changeVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename CHVAR>
void VariableChange<MODEL, CHVAR>::changeVarInverse(const State_ & x1, State_ & x2) const {
  Log::trace() << "VariableChange<MODEL, CHVAR>::changeVarInverse starting" << std::endl;
  util::Timer timer(classname(), "changeVarInverse");
  chvar_->changeVarInverse(x1.state(), x2.state());
  Log::trace() << "VariableChange<MODEL, CHVAR>::changeVarInverse done" << std::endl;
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
