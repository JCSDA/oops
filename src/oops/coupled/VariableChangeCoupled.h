/*
 * (C) Copyright 2022-2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/coupled/UtilsCoupled.h"
#include "oops/interface/VariableChange.h"
#include "oops/util/Printable.h"

namespace oops {
  class Variables;

// -----------------------------------------------------------------------------
/// Change of variables for coupled state
template <typename MODEL1, typename MODEL2>
class VariableChangeCoupled : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>  Geometry_;
  typedef StateCoupled<MODEL1, MODEL2>     State_;
  typedef VariableChange<MODEL1>           VariableChange1_;
  typedef VariableChange<MODEL2>           VariableChange2_;

 public:
  VariableChangeCoupled(const eckit::Configuration &, const Geometry_ &);

  /// Perform transforms
  void changeVar(State_ &, const oops::Variables &) const;
  void changeVarInverse(State_ &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<VariableChange1_> varchg1_;
  std::unique_ptr<VariableChange2_> varchg2_;
  const std::vector<Variables> availableVars_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
VariableChangeCoupled<MODEL1, MODEL2>::VariableChangeCoupled(const eckit::Configuration & config,
                                                             const Geometry_ & geometry)
  : availableVars_(geometry.variables())
{
  if (geometry.isParallel()) {
    throw eckit::NotImplemented(Here());
  }
  const eckit::LocalConfiguration conf1 = config.getSubConfiguration(MODEL1::name());;
  const eckit::LocalConfiguration conf2 = config.getSubConfiguration(MODEL2::name());;
  varchg1_ = std::make_unique<VariableChange1_>(conf1, geometry.geometry1());
  varchg2_ = std::make_unique<VariableChange2_>(conf2, geometry.geometry2());
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void VariableChangeCoupled<MODEL1, MODEL2>::changeVar(State_ & xx,
                                            const oops::Variables & vars) const {
  // decide what variables are provided by what model
  std::vector<Variables> splitvars = splitVariables(vars, availableVars_);
  varchg1_->changeVar(xx.state1(), splitvars[0]);
  varchg2_->changeVar(xx.state2(), splitvars[1]);
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void VariableChangeCoupled<MODEL1, MODEL2>::changeVarInverse(State_ & xx,
                                            const oops::Variables & vars) const {
  // decide what variables are provided by what model
  std::vector<Variables> splitvars = splitVariables(vars, availableVars_);
  varchg1_->changeVarInverse(xx.state1(), splitvars[0]);
  varchg2_->changeVarInverse(xx.state2(), splitvars[1]);
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void VariableChangeCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  os << "Coupled VariableChange (" << MODEL1::name() << " and " <<
        MODEL2::name() << ")" << std::endl;
  os << *varchg1_ << std::endl << *varchg2_;
}

// -----------------------------------------------------------------------------

}  // namespace oops
