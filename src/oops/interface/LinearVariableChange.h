/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_LINEARVARIABLECHANGE_H_
#define OOPS_INTERFACE_LINEARVARIABLECHANGE_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearVariableChangeParametersBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/GenericParameters.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief MODEL-agnostic part of the linear variable change
///
/// Note: each implementation should typedef `LinVariableChangeParams_` to the name of a
/// subclass of oops::LinearVariableChangeParametersBase holding its configuration settings and
/// provide a constructor with the following signature:
///
///     LinearVariableChange(const Geometry_ &, const LinVariableChangeParams_ &);
template <typename MODEL>
class LinearVariableChange {
  typedef typename MODEL::LinearVariableChange LinearVariableChange_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  /// A subclass of oops::LinearVariableChangeParametersBase holding the configuration settings of
  /// the variable change.
  /// Defined as LinearVariableChange_::Parameters_ if LinearVariableChange_ defines a Parameters_
  /// type; otherwise as GenericLinearVariableChangeParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<LinearVariableChange_,
                                    GenericLinearVariableChangeParameters> Parameters_;
  static const std::string classname() {return "oops::LinearVariableChange";}

  LinearVariableChange(const Geometry_ &, const Parameters_ &);
  LinearVariableChange(const Geometry_ &, const eckit::Configuration &);

  virtual ~LinearVariableChange();

  void changeVarTraj(const State_ &, const Variables &);
  void changeVarTL(Increment_ &, const Variables &) const;
  void changeVarInverseTL(Increment_ &, const Variables &) const;
  void changeVarAD(Increment_ &, const Variables &) const;
  void changeVarInverseAD(Increment_ &, const Variables &) const;

 private:
  void print(std::ostream &) const;

  std::unique_ptr<LinearVariableChange_> chvar_;
};

// =============================================================================

template<typename MODEL>
LinearVariableChange<MODEL>::LinearVariableChange(const Geometry_ & resol,
                                                  const eckit::Configuration & conf) : chvar_() {
  Log::trace() << "LinearVariableChange<MODEL>::LinearVariableChange starting" << std::endl;
  util::Timer timer(classname(), "LinearVariableChange");
  chvar_.reset(new LinearVariableChange_(resol.geometry(), conf));
  Log::trace() << "LinearVariableChange<MODEL>::LinearVariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearVariableChange<MODEL>::LinearVariableChange(const Geometry_ & resol,
                                                  const Parameters_ & parameters)
  :  LinearVariableChange<MODEL>::LinearVariableChange(resol, parameters.toConfiguration())
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearVariableChange<MODEL>::~LinearVariableChange() {
  Log::trace() << "LinearVariableChange<MODEL>::~LinearVariableChange starting" << std::endl;
  util::Timer timer(classname(), "~LinearVariableChange");
  chvar_.reset();
  Log::trace() << "LinearVariableChange<MODEL>::~LinearVariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::changeVarTL(Increment_ & dx, const Variables & vars) const {
  Log::trace() << "LinearVariableChange<MODEL>::changeVarTL starting" << std::endl;
  util::Timer timer(classname(), "changeVarTL");
  chvar_->changeVarTL(dx.increment(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::changeVarTL done" << std::endl;
}
// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::changeVarInverseTL(Increment_ & dx,
                                                     const Variables & vars) const {
  Log::trace() << "LinearVariableChange<MODEL>::changeVarInverseTL starting" << std::endl;
  util::Timer timer(classname(), "changeVarInverseTL");
  chvar_->changeVarInverseTL(dx.increment(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::changeVarInverseTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::changeVarAD(Increment_ & dx, const Variables & vars) const {
  Log::trace() << "LinearVariableChange<MODEL>::changeVarAD starting" << std::endl;
  util::Timer timer(classname(), "changeVarAD");
  chvar_->changeVarAD(dx.increment(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::changeVarAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::changeVarInverseAD(Increment_ & dx,
                                                     const Variables & vars) const {
  Log::trace() << "LinearVariableChange<MODEL>::changeVarInverseAD starting" << std::endl;
  util::Timer timer(classname(), "changeVarInverseAD");
  chvar_->changeVarInverseAD(dx.increment(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::changeVarInverseAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::changeVarTraj(const State_ & xFirstGuess,
                                                const Variables & vars) {
  Log::trace() << "LinearVariableChange<MODEL>::changeVarTraj starting" << std::endl;
  util::Timer timer(classname(), "changeVarTraj");
  chvar_->changeVarTraj(xFirstGuess.state(), vars);
  Log::trace() << "LinearVariableChange<MODEL>::changeVarTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearVariableChange<MODEL>::print(std::ostream & os) const {
  Log::trace() << "LinearVariableChange<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *chvar_;
  Log::trace() << "LinearVariableChange<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEARVARIABLECHANGE_H_
