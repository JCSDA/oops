/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2021 UCAR.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTFCT4DVAR_H_
#define OOPS_ASSIMILATION_COSTFCT4DVAR_H_

#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJcDFI.h"
#include "oops/assimilation/CostJo.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearModel.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"

namespace oops {

/// Strong Constraint 4D-Var Cost Function
/*!
 * This class is not really necessary since it is only a special
 * case of the more general weak constraint 4D-Var cost function
 * with one sub-window. It is provided for readability.
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class CostFct4DVar : public CostFunction<MODEL, OBS> {
  typedef Increment<MODEL>                Increment_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef ControlVariable<MODEL, OBS>     CtrlVar_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;
  typedef Model<MODEL>                    Model_;
  typedef LinearModel<MODEL>              LinearModel_;

 public:
  CostFct4DVar(const eckit::Configuration &, const eckit::mpi::Comm &);
  ~CostFct4DVar() {}

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>,
              const bool idModel = false) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>,
              const bool idModel = false) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

 protected:
  const Geometry_ & geometry() const override {return resol_;}

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const override;

  CostJb3D<MODEL, OBS> * newJb(const eckit::Configuration &, const Geometry_ &) const override;
  CostJo<MODEL, OBS>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &, CtrlVar_ &, CtrlVar_ &,
                   PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) override;

  const eckit::mpi::Comm & comm_;
  const util::TimeWindow timeWindow_;
  const Geometry_ resol_;
  Model_ model_;
  const Variables ctlvars_;
  std::shared_ptr<LinearModel_> tlm_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFct4DVar<MODEL, OBS>::CostFct4DVar(const eckit::Configuration & config,
                                       const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(), comm_(comm),
    timeWindow_(config.getSubConfiguration("time window")),
    resol_(eckit::LocalConfiguration(config, "geometry"), comm),
    model_(resol_, eckit::LocalConfiguration(config, "model")),
    ctlvars_(config.getStringVector("analysis variables")), tlm_()
{
  Log::trace() << "CostFct4DVar:CostFct4DVar start" << std::endl;

  this->setupTerms(config);
  // ASSERT(ctlvars_ <= this->jb().getBackground().state().variables());
  Log::info() << "4DVar window: " << timeWindow_ << std::endl;
  Log::trace() << "CostFct4DVar::CostFct4DVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb3D<MODEL, OBS> * CostFct4DVar<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                  const Geometry_ & resol) const {
  Log::trace() << "CostFct4DVar::newJb start" << std::endl;
  return new CostJb3D<MODEL, OBS>(timeWindow_.start(), jbConf, resol, ctlvars_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFct4DVar<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFct4DVar::newJo start" << std::endl;
  return new CostJo<MODEL, OBS>(joConf, comm_, timeWindow_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFct4DVar<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                           const Geometry_ & resol) const {
  Log::trace() << "CostFct4DVar::newJc start" << std::endl;
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(timeWindow_.midpoint());
  return new CostJcDFI<MODEL, OBS>(jcdfi, resol, vt, timeWindow_.length());
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFct4DVar::runNL start" << std::endl;
  ASSERT(xx.states().is_3d());
  ASSERT(xx.state().validTime() == timeWindow_.start());
  model_.forecast(xx.state(), xx.modVar(), timeWindow_.length(), post);
  ASSERT(xx.state().validTime() == timeWindow_.end());
  Log::trace() << "CostFct4DVar::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::doLinearize(const Geometry_ & resol,
                                           const eckit::Configuration & innerConf,
                                           CtrlVar_ & bg, CtrlVar_ & fg,
                                           PostProcessor<State_> & pp,
                                           PostProcessorTLAD<MODEL> & pptraj) {
  Log::trace() << "CostFct4DVar::doLinearize start" << std::endl;
  ASSERT(bg.states().is_3d());
  ASSERT(fg.states().is_3d());
  eckit::LocalConfiguration lmConf(innerConf, "linear model");
// Setup linear model (and trajectory)
  tlm_.reset(new LinearModel_(resol, lmConf));
  pp.enrollProcessor(new TrajectorySaver<MODEL>(lmConf, resol, fg.modVar(), tlm_, pptraj));

  Log::trace() << "CostFct4DVar::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                      PostProcessorTLAD<MODEL> & cost,
                                      PostProcessor<Increment_> post,
                                      const bool idModel) const {
  Log::trace() << "CostFct4DVar::runTLM start" << std::endl;
  ASSERT(dx.states().is_3d());
  ASSERT(dx.state().validTime() == timeWindow_.start());

  tlm_->forecastTL(dx.state(), dx.modVar(), timeWindow_.length(), post, cost, idModel);
  ASSERT(dx.state().validTime() == timeWindow_.end());
  Log::trace() << "CostFct4DVar::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFct4DVar::zeroAD start" << std::endl;
  ASSERT(dx.states().is_3d());
  dx.state().zero(timeWindow_.end());
  dx.modVar().zero();
  dx.obsVar().zero();
  Log::trace() << "CostFct4DVar::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                      PostProcessorTLAD<MODEL> & cost,
                                      PostProcessor<Increment_> post,
                                      const bool idModel) const {
  Log::trace() << "CostFct4DVar::runADJ start" << std::endl;
  ASSERT(dx.states().is_3d());
  ASSERT(dx.state().validTime() == timeWindow_.end());

  tlm_->forecastAD(dx.state(), dx.modVar(), timeWindow_.length(), post, cost, idModel);

  ASSERT(dx.state().validTime() == timeWindow_.start());
  Log::trace() << "CostFct4DVar::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                       PostProcessor<Increment_> &) const {
  Log::trace() << "CostFct4DVar::addIncr start" << std::endl;
  ASSERT(xx.states().is_3d());
  ASSERT(dx.states().is_3d());
  xx.state() += dx.state();
  Log::trace() << "CostFct4DVar::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT4DVAR_H_
