/*
 * (C) Copyright 2022- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/interface/VariableChange.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace oops {
  class Variables;

/// Parameters describing a variable change for coupled state
template <typename MODEL1, typename MODEL2>
class VariableChangeCoupledParameters : public VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeCoupledParameters,
                           VariableChangeParametersBase)
  typedef typename VariableChange<MODEL1>::Parameters_ Parameters1_;
  typedef typename VariableChange<MODEL2>::Parameters_ Parameters2_;
 public:
  RequiredParameter<Parameters1_> varchg1{MODEL1::name().c_str(), this};
  RequiredParameter<Parameters2_> varchg2{MODEL2::name().c_str(), this};
  RequiredParameter<Variables> vars1{std::string(MODEL1::name() + " variables").c_str(),
          "variables that the first model should provide, have to be different "
          "from the variables that the second model provides", this};
  RequiredParameter<Variables> vars2{std::string(MODEL2::name() + " variables").c_str(),
          "variables that the second model should provide, have to be different "
          "from the variables that the first model provides", this};
};

// -----------------------------------------------------------------------------
/// Change of variables for coupled state
template <typename MODEL1, typename MODEL2>
class VariableChangeCoupled : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>  Geometry_;
  typedef StateCoupled<MODEL1, MODEL2>     State_;
  typedef VariableChange<MODEL1>           VariableChange1_;
  typedef VariableChange<MODEL2>           VariableChange2_;

 public:
  typedef VariableChangeCoupledParameters<MODEL1, MODEL2> Parameters_;

  VariableChangeCoupled(const Parameters_ &, const Geometry_ &);

  /// Perform transforms
  void changeVar(State_ &, const oops::Variables &) const;
  void changeVarInverse(State_ &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<VariableChange1_> varchg1_;
  std::unique_ptr<VariableChange2_> varchg2_;
  const Variables vars1_;  ///< variables that model1 should provide
  const Variables vars2_;  ///< variables that model2 should provide
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
VariableChangeCoupled<MODEL1, MODEL2>::VariableChangeCoupled(
      const Parameters_ & params, const Geometry_ & geometry)
  : vars1_(params.vars1), vars2_(params.vars2) {
  if (geometry.isParallel()) {
    throw eckit::NotImplemented(Here());
  }
  // check that the same variable isn't specified in both models'
  // variables
  Variables commonvars = vars1_;
  commonvars.intersection(vars2_);
  if (commonvars.size() > 0) {
    throw eckit::BadParameter("Variables for different components of coupled "
          "variable change can not overlap", Here());
  }

  varchg1_ = std::make_unique<VariableChange1_>(params.varchg1, geometry.geometry1());
  varchg2_ = std::make_unique<VariableChange2_>(params.varchg2, geometry.geometry2());
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void VariableChangeCoupled<MODEL1, MODEL2>::changeVar(State_ & xx,
                                            const oops::Variables & vars) const {
  // decide what variables are provided by what model
  Variables outvars1 = vars;
  outvars1.intersection(vars1_);
  Variables outvars2 = vars;
  outvars2.intersection(vars2_);
  // check that all variables are accounted for
  Variables alloutvars = outvars1;
  alloutvars += outvars2;
  if (alloutvars != vars) {
    throw eckit::UserError("Not all variables can be provided by the coupled variable change");
  }
  varchg1_->changeVar(xx.state1(), outvars1);
  varchg2_->changeVar(xx.state2(), outvars2);
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void VariableChangeCoupled<MODEL1, MODEL2>::changeVarInverse(State_ & xx,
                                            const oops::Variables & vars) const {
  Variables outvars1 = vars;
  outvars1.intersection(vars1_);
  Variables outvars2 = vars;
  outvars2.intersection(vars2_);
  // check that all variables are accounted for
  Variables alloutvars = outvars1;
  alloutvars += outvars2;
  if (alloutvars != vars) {
    throw eckit::UserError("Not all variables can be provided by the coupled variable change");
  }
  varchg1_->changeVarInverse(xx.state1(), outvars1);
  varchg2_->changeVarInverse(xx.state2(), outvars2);
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void VariableChangeCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  os << "Coupled VariableChange (" << MODEL1::name() << " and " <<
        MODEL2::name() << ")" << std::endl;
  os << *varchg1_ << std::endl << *varchg2_;
}

}  // namespace oops
