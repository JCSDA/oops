/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "oops/base/LinearModelBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"

namespace oops {
namespace interface {

// -----------------------------------------------------------------------------

template<typename MODEL>
class LinearModel : public oops::LinearModelBase<MODEL> {
  typedef typename MODEL::LinearModel     LinearModel_;
  typedef oops::Geometry<MODEL>           Geometry_;
  typedef oops::Increment<MODEL>          Increment_;
  typedef oops::ModelAuxControl<MODEL>    ModelAuxCtl_;
  typedef oops::ModelAuxIncrement<MODEL>  ModelAuxInc_;
  typedef oops::State<MODEL>              State_;

 public:
  static const std::string classname() {return "oops::interface::LinearModel";}

  LinearModel(const Geometry_ &, const eckit::Configuration &);
  ~LinearModel();

// Set the linearization trajectory
  void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &) override;

// TL forecast
  void initializeTL(Increment_ &) const override;
  void stepTL(Increment_ &, const ModelAuxInc_ &) const override;
  void finalizeTL(Increment_ &) const override;

// AD forecast
  void initializeAD(Increment_ &) const override;
  void stepAD(Increment_ &, ModelAuxInc_ &) const override;
  void finalizeAD(Increment_ &) const override;

// Information and diagnostics
  const util::Duration & timeResolution() const override {return tlm_->timeResolution();}
  const util::Duration & stepTrajectory() const override {return tlm_->stepTrajectory();}

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<LinearModel_> tlm_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearModel<MODEL>::LinearModel(const Geometry_ & resol, const eckit::Configuration & conf)
  : tlm_()
{
  Log::trace() << "LinearModel<MODEL>::LinearModel starting" << std::endl;
  util::Timer timer(classname(), "LinearModel");
  tlm_.reset(new LinearModel_(resol.geometry(), conf));
  Log::trace() << "LinearModel<MODEL>::LinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearModel<MODEL>::~LinearModel() {
  Log::trace() << "LinearModel<MODEL>::~LinearModel starting" << std::endl;
  util::Timer timer(classname(), "~LinearModel");
  tlm_.reset();
  Log::trace() << "LinearModel<MODEL>::~LinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::setTrajectory(const State_ & xx, State_ & xlr, const ModelAuxCtl_ & maux) {
  Log::trace() << "LinearModel<MODEL>::setTrajectory starting" << std::endl;
  util::Timer timer(classname(), "setTrajectory");
  tlm_->setTrajectory(xx.state(), xlr.state(), maux.modelauxcontrol());
  Log::trace() << "LinearModel<MODEL>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::initializeTL(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::initializeTL starting" << std::endl;
  util::Timer timer(classname(), "initializeTL");
  tlm_->initializeTL(dx.increment());
  Log::trace() << "LinearModel<MODEL>::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::stepTL(Increment_ & dx, const ModelAuxInc_ & merr) const {
  Log::trace() << "LinearModel<MODEL>::stepTL starting" << std::endl;
  util::Timer timer(classname(), "stepTL");
  tlm_->stepTL(dx.increment(), merr.modelauxincrement());
  Log::trace() << "LinearModel<MODEL>::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::finalizeTL starting" << std::endl;
  util::Timer timer(classname(), "finalizeTL");
  tlm_->finalizeTL(dx.increment());
  Log::trace() << "LinearModel<MODEL>::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
void LinearModel<MODEL>::initializeAD(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::initializeAD starting" << std::endl;
  util::Timer timer(classname(), "initializeAD");
  tlm_->initializeAD(dx.increment());
  Log::trace() << "LinearModel<MODEL>::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::stepAD(Increment_ & dx, ModelAuxInc_ & merr) const {
  Log::trace() << "LinearModel<MODEL>::stepAD starting" << std::endl;
  util::Timer timer(classname(), "stepAD");
  tlm_->stepAD(dx.increment(), merr.modelauxincrement());
  Log::trace() << "LinearModel<MODEL>::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::finalizeAD starting" << std::endl;
  util::Timer timer(classname(), "finalizeAD");
  tlm_->finalizeAD(dx.increment());
  Log::trace() << "LinearModel<MODEL>::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::print(std::ostream & os) const {
  Log::trace() << "LinearModel<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *tlm_;
  Log::trace() << "LinearModel<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace interface
}  // namespace oops
