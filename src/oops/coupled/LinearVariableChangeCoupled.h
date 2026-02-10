/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/IncrementCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Coupled linear variable change assuming the linear variable changes
///        for the two components are not interacting with each other.
///
template <typename MODEL1, typename MODEL2>
class LinearVariableChangeCoupled {
  typedef GeometryCoupled<MODEL1, MODEL2>  GeometryCoupled_;
  typedef IncrementCoupled<MODEL1, MODEL2> IncrementCoupled_;
  typedef StateCoupled<MODEL1, MODEL2>     StateCoupled_;

 public:
  static const std::string classname() {return "oops::LinearVariableChangeCoupled";}

  LinearVariableChangeCoupled(const GeometryCoupled_ &, const eckit::Configuration &);

  virtual ~LinearVariableChangeCoupled() = default;

  void changeVarTraj(const StateCoupled_ &, const Variables &);
  void changeVarTL(IncrementCoupled_ &, const Variables &) const;
  void changeVarInverseTL(IncrementCoupled_ &, const Variables &) const;
  void changeVarAD(IncrementCoupled_ &, const Variables &) const;
  void changeVarInverseAD(IncrementCoupled_ &, const Variables &) const;

 private:
  std::unique_ptr<LinearVariableChange<MODEL1>> chvar1_;
  std::unique_ptr<LinearVariableChange<MODEL2>> chvar2_;
  const std::vector<Variables> availableVars_;
};

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
LinearVariableChangeCoupled<MODEL1, MODEL2>::LinearVariableChangeCoupled(
    const GeometryCoupled_ & resol,
    const eckit::Configuration & conf)
  : chvar1_(), chvar2_(), availableVars_(resol.variables()) {
  Log::trace() << "LinearVariableChangeCoupled::LinearVariableChangeCoupled starting" << std::endl;
  util::Timer timer(classname(), "LinearVariableChangeCoupled");
  const eckit::LocalConfiguration conf1 = conf.getSubConfiguration(MODEL1::name());;
  const eckit::LocalConfiguration conf2 = conf.getSubConfiguration(MODEL2::name());;
  chvar1_.reset(new LinearVariableChange<MODEL1>(resol.geometry1(), conf1));
  chvar2_.reset(new LinearVariableChange<MODEL2>(resol.geometry2(), conf2));
  Log::trace() << "LinearVariableChangeCoupled::LinearVariableChangeCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void LinearVariableChangeCoupled<MODEL1, MODEL2>::changeVarTraj(
    const StateCoupled_ & xFirstGuess, const Variables & vars) {
  Log::trace() << "LinearVariableChangeCoupled::changeVarTraj starting" << std::endl;
  util::Timer timer(classname(), "changeVarTraj");
  std::vector<Variables> splitvars = splitVariables(vars, availableVars_);
  chvar1_->changeVarTraj(xFirstGuess.state1(), splitvars[0]);
  chvar2_->changeVarTraj(xFirstGuess.state2(), splitvars[1]);
  Log::trace() << "LinearVariableChangeCoupled::changeVarTraj done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void LinearVariableChangeCoupled<MODEL1, MODEL2>::changeVarTL(
    IncrementCoupled_ & dx, const Variables & vars) const {
  Log::trace() << "LinearVariableChangeCoupled::changeVarTL starting" << std::endl;
  util::Timer timer(classname(), "changeVarTL");
  std::vector<Variables> splitvars = splitVariables(vars, availableVars_);
  chvar1_->changeVarTL(dx.increment1(), splitvars[0]);
  chvar2_->changeVarTL(dx.increment2(), splitvars[1]);
  Log::trace() << "LinearVariableChangeCoupled::changeVarTL done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void LinearVariableChangeCoupled<MODEL1, MODEL2>::changeVarInverseTL(
    IncrementCoupled_ & dx, const Variables & vars) const {
  Log::trace() << "LinearVariableChangeCoupled::changeVarInverseTL starting" << std::endl;
  util::Timer timer(classname(), "changeVarInverseTL");
  std::vector<Variables> splitvars = splitVariables(vars, availableVars_);
  chvar1_->changeVarInverseTL(dx.increment1(), splitvars[0]);
  chvar2_->changeVarInverseTL(dx.increment2(), splitvars[1]);
  Log::trace() << "LinearVariableChangeCoupled::changeVarInverseTL done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void LinearVariableChangeCoupled<MODEL1, MODEL2>::changeVarAD(
    IncrementCoupled_ & dx, const Variables & vars) const {
  Log::trace() << "LinearVariableChangeCoupled::changeVarAD starting" << std::endl;
  util::Timer timer(classname(), "changeVarAD");
  std::vector<Variables> splitvars = splitVariables(vars, availableVars_);
  chvar1_->changeVarAD(dx.increment1(), splitvars[0]);
  chvar2_->changeVarAD(dx.increment2(), splitvars[1]);
  Log::trace() << "LinearVariableChangeCoupled::changeVarAD done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void LinearVariableChangeCoupled<MODEL1, MODEL2>::changeVarInverseAD(
    IncrementCoupled_ & dx, const Variables & vars) const {
  Log::trace() << "LinearVariableChangeCoupled::changeVarInverseAD starting" << std::endl;
  util::Timer timer(classname(), "changeVarInverseAD");
  std::vector<Variables> splitvars = splitVariables(vars, availableVars_);
  chvar1_->changeVarInverseAD(dx.increment1(), splitvars[0]);
  chvar2_->changeVarInverseAD(dx.increment2(), splitvars[1]);
  Log::trace() << "LinearVariableChangeCoupled::changeVarInverseAD done" << std::endl;
}

}  // namespace oops
