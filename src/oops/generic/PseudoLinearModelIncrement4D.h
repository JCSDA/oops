/*
 * (C) Copyright 2023- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/generic/LinearModelBase.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic implementation of the pseudo linear model initialized with 4D Increment
/// (steps through time by stepping through states in 4D Increment)
template <typename MODEL>
class PseudoLinearModelIncrement4D : public LinearModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef Increment<MODEL>         Increment_;
  typedef Increment4D<MODEL>       Increment4D_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef ModelAuxIncrement<MODEL> ModelAuxInc_;
  typedef State<MODEL>             State_;
  typedef State4D<MODEL>           State4D_;

 public:
  static const std::string classname() {return "oops::PseudoLinearModelIncrement4D";}

  /// Initialize pseudo model with \p inc - 4D increment to loop through
  /// in the model run and \p tstep - time resolution of the model
  PseudoLinearModelIncrement4D(const Increment4D_ & inc,
                    const util::Duration & tstep = util::Duration(0));

/// initialize tangent linear forecast
  void initializeTL(Increment_ &) const override;
/// one tangent linear forecast step
  void stepTL(Increment_ &, const ModelAuxInc_ &) const override;
/// finalize tangent linear forecast
  void finalizeTL(Increment_ &) const override;

/// initialize adjoint forecast
  void initializeAD(Increment_ &) const override;
/// one adjoint forecast step
  void stepAD(Increment_ &, ModelAuxInc_ &) const override;
/// finalize adjoint forecast
  void finalizeAD(Increment_ &) const override;

/// set trajectory
  void setTrajectory(const State_ &, State_ &, const ModelAux_ &) override {}
  // void setTrajectory(const State_ &, State_ &, const ModelAux_ &);

  /// linear model time step
  const util::Duration & timeResolution() const override {return tstep_;}
  /// linear model variables
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const override;

  /// Reference to 4D state that is used in the model
  const Increment4D_ & inc4d_;
  /// Model's time resolution
  util::Duration   tstep_;
  /// Index of the current increment
  mutable size_t currentinc_;
  oops::Variables vars_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
PseudoLinearModelIncrement4D<MODEL>::PseudoLinearModelIncrement4D(const Increment4D_ & inc4d,
                                              const util::Duration & tstep)
  : inc4d_(inc4d), tstep_(tstep) {
  const std::vector<util::DateTime> validTimes = inc4d_.validTimes();
  vars_ = inc4d_.variables();
  if (validTimes.size() > 1) tstep_ = validTimes[1] - validTimes[0];
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>::PseudoLinearModelIncrement4D done"
               << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoLinearModelIncrement4D<MODEL>::initializeTL(Increment_ & dx) const {
  currentinc_ = 0;
  dx = inc4d_[currentinc_];
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoLinearModelIncrement4D<MODEL>::stepTL(Increment_ & dx,
                                                 const ModelAuxInc_ & mainc) const {
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>:step Starting " << std::endl;
  dx = inc4d_[currentinc_++];
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoLinearModelIncrement4D<MODEL>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoLinearModelIncrement4D<MODEL>::initializeAD(Increment_ & dx) const {
  currentinc_ = inc4d_.size()-1;
  dx = inc4d_[currentinc_];
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoLinearModelIncrement4D<MODEL>::stepAD(Increment_ & dx, ModelAuxInc_ & merr) const {
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>:step Starting " << std::endl;
  dx = inc4d_[currentinc_--];
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoLinearModelIncrement4D<MODEL>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "PseudoLinearModelIncrement4D<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoLinearModelIncrement4D<MODEL>::print(std::ostream & os) const {
  os << "Pseudo model stepping through 4D state with " << tstep_ << " time resolution";
}

}  // namespace oops
