/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_VARIABLECHANGE_H_
#define OOPS_INTERFACE_VARIABLECHANGE_H_

#include <memory>
#include <string>

#include "oops/base/VariableChangeBase.h"
#include "oops/base/VariableChangeParametersBase.h"
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

/// \brief Encapsulates the nonlinear variable change
/// Note: to see methods that need to be implemented in the implementation,
/// see VariableChangeBase class.
template <typename MODEL>
class VariableChange : public util::Printable,
                       private util::ObjectCounter<VariableChange<MODEL> >  {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef GenericVariableChangeBase<MODEL>  VariableChangeBase_;

 public:
  static const std::string classname() {return "oops::VariableChange";}

  VariableChange(const Geometry_ &, const VariableChangeParametersBase &);
  VariableChange(const Geometry_ &, const eckit::Configuration &);
  virtual ~VariableChange();
  VariableChange(const VariableChange &) = delete;
  VariableChange(VariableChange &&) = default;
  const VariableChange & operator=(const VariableChange &) = delete;
  VariableChange & operator=(VariableChange &&) = default;

  /// set input variables for variable transform
  void setInputVariables(const Variables & vars) { varin_.reset(new Variables(vars)); }
  /// set output variables for variable transform
  void setOutputVariables(const Variables & vars) { varout_.reset(new Variables(vars)); }

  /// change variable from state \p xin to \p xout
  void changeVar(const State_ & xin, State_ & xout) const;
  /// inverse of changeVar, change variables back from \p xout to \p xin
  void changeVarInverse(const State_ & xout, State_ & xin) const;

  /// return change of variable \p xin
  State_ changeVar(const State_ & xin) const;
  /// return inverse of variable change applied to \p xout
  State_ changeVarInverse(const State_ & xout) const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<VariableChangeBase_> chvar_;  /// pointer to the VariableChange implementation
  std::unique_ptr<Variables> varin_;            /// input variables
  std::unique_ptr<Variables> varout_;           /// output variables
};

// =============================================================================

template<typename MODEL>
VariableChange<MODEL>::VariableChange(const Geometry_ & geom,
                                      const VariableChangeParametersBase & params)
  : chvar_()
{
  Log::trace() << "VariableChange<MODEL>::VariableChange starting" << std::endl;
  util::Timer timer(classname(), "VariableChange");
  if (params.inputVariables.value() != boost::none) {
    varin_.reset(new Variables(*params.inputVariables.value()));
    Log::trace() << "VariableChange<MODEL>::VariableChange input variables: "
                 << *varin_ << std::endl;
  }
  if (params.outputVariables.value() != boost::none) {
    varout_.reset(new Variables(*params.outputVariables.value()));
    Log::trace() << "VariableChange<MODEL>::VariableChange output variables: "
                 << *varout_ << std::endl;
  }
  chvar_.reset(VariableChangeFactory<MODEL>::create(geom, params));
  Log::trace() << "VariableChange<MODEL>::VariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
VariableChange<MODEL>::VariableChange(const Geometry_ & geom, const eckit::Configuration & conf)
  : VariableChange(geom,
    validateAndDeserialize<VariableChangeParametersWrapper<MODEL>>(conf).variableChangeParameters)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
VariableChange<MODEL>::~VariableChange() {
  Log::trace() << "VariableChange<MODEL>::~VariableChange starting" << std::endl;
  util::Timer timer(classname(), "~VariableChange");
  chvar_.reset();
  Log::trace() << "VariableChange<MODEL>::~VariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void VariableChange<MODEL>::changeVar(const State_ & x1, State_ & x2) const {
  Log::trace() << "VariableChange<MODEL>::changeVar starting" << std::endl;
  util::Timer timer(classname(), "changeVar");
  chvar_->changeVar(x1, x2);
  Log::trace() << "VariableChange<MODEL>::changeVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void VariableChange<MODEL>::changeVarInverse(const State_ & x1, State_ & x2) const {
  Log::trace() << "VariableChange<MODEL>::changeVarInverse starting" << std::endl;
  util::Timer timer(classname(), "changeVarInverse");
  chvar_->changeVarInverse(x1, x2);
  Log::trace() << "VariableChange<MODEL>::changeVarInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> VariableChange<MODEL>::changeVar(const State_ & xin) const {
  Log::trace() << "VariableChange<MODEL>::changeVar starting" << std::endl;
  ASSERT(varout_);
  State_ xout(xin.geometry(), *varout_, xin.validTime());
  this->changeVar(xin, xout);
  Log::trace() << "VariableChange<MODEL>::changeVar done" << std::endl;
  return xout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> VariableChange<MODEL>::changeVarInverse(const State_ & xin) const {
  Log::trace() << "VariableChange<MODEL>::changeVarInverse starting" << std::endl;
  ASSERT(varin_);
  State_ xout(xin.geometry(), *varin_, xin.validTime());
  this->changeVarInverse(xin, xout);
  Log::trace() << "VariableChange<MODEL>::changeVarInverse done" << std::endl;
  return xout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void VariableChange<MODEL>::print(std::ostream & os) const {
  Log::trace() << "VariableChange<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *chvar_;
  if (varin_)  os << std::endl << "Variable change from: " << *varin_;
  if (varout_) os << std::endl << "Variable change to: " << *varout_;
  if (varin_ || varout_) os << std::endl;
  Log::trace() << "VariableChange<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_VARIABLECHANGE_H_
