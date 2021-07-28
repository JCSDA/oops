/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_LINEARMODELID_H_
#define OOPS_GENERIC_LINEARMODELID_H_

#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearModelBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Encapsulates the linear forecast model.
/*!
 * Generic implementation of the identity linear model.
 */

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearModelId : public LinearModelBase<MODEL> {
  typedef typename MODEL::Increment         Increment_;
  typedef typename MODEL::Geometry          Geometry_;
  typedef typename MODEL::ModelAuxControl   ModelAux_;
  typedef typename MODEL::ModelAuxIncrement ModelAuxIncr_;
  typedef typename MODEL::State             State_;

 public:
  static const std::string classname() {return "oops::LinearModelId";}

  LinearModelId(const Geometry_ &, const eckit::Configuration &);
  ~LinearModelId();

// Set the linearization trajectory
  void setTrajectory(const State_ &, State_ &, const ModelAux_ &) override;

// Run the TL forecast
  void initializeTL(Increment_ &) const override;
  void stepTL(Increment_ & dx, const ModelAuxIncr_ &) const override;
  void finalizeTL(Increment_ &) const override;

// Run the AD forecast:
  void initializeAD(Increment_ &) const override;
  void stepAD(Increment_ & dx, ModelAuxIncr_ &) const override;
  void finalizeAD(Increment_ &) const override;

// Information and diagnostics
  const util::Duration & timeResolution() const override {return tstep_;}
  const oops::Variables & variables() const override {return vars_;}
  void print(std::ostream &) const override {}

 private:
  const util::Duration tstep_;
  const Variables vars_;
};

// =============================================================================

template<typename MODEL>
LinearModelId<MODEL>::LinearModelId(const Geometry_ &, const eckit::Configuration & tlConf)
  : tstep_(util::Duration(tlConf.getString("tstep"))), vars_()
{
  Log::trace() << "LinearModelId<MODEL>::LinearModelId" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearModelId<MODEL>::~LinearModelId() {
  Log::trace() << "LinearModelId<MODEL>::~LinearModelId" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelId<MODEL>::setTrajectory(const State_ &, State_ &, const ModelAux_ &) {
  Log::trace() << "LinearModelId<MODEL>::setTrajectory" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelId<MODEL>::initializeTL(Increment_ &) const {
  Log::trace() << "LinearModelId<MODEL>::initializeTL" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelId<MODEL>::stepTL(Increment_ & dx, const ModelAuxIncr_ &) const {
  dx.updateTime(tstep_);
  Log::trace() << "LinearModelId<MODEL>::stepTL" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelId<MODEL>::finalizeTL(Increment_ &) const {
  Log::trace() << "LinearModelId<MODEL>::finalizeTL" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelId<MODEL>::initializeAD(Increment_ &) const {
  Log::trace() << "LinearModelId<MODEL>::initializeAD" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelId<MODEL>::stepAD(Increment_ & dx, ModelAuxIncr_ &) const {
  dx.updateTime(-tstep_);
  Log::trace() << "LinearModelId<MODEL>::stepAD" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelId<MODEL>::finalizeAD(Increment_ &) const {
  Log::trace() << "LinearModelId<MODEL>::finalizeAD" << std::endl;
}

// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_GENERIC_LINEARMODELID_H_
