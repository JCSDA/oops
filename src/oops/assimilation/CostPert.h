/*
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_COSTPERT_H_
#define OOPS_ASSIMILATION_COSTPERT_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFct3DVar.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/StateSaver.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

/// 3D-Var Cost Function for the Pert members of Control-Pert EDA

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS, class J> class CostPert : public J
{
    typedef Increment<MODEL>                Increment_;
    typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
    typedef ControlVariable<MODEL, OBS>     CtrlVar_;
    typedef CostFct3DVar<MODEL, OBS>        CostFct3DVar_;
    typedef CostFunction<MODEL, OBS>        CostFct_;
    typedef CostTermBase<MODEL, OBS>        CostBase_;
    typedef Geometry<MODEL>                 Geometry_;
    typedef ModelAuxControl<MODEL>          ModelAuxControl_;
    typedef ObsAuxControls<OBS>             ObsAuxControls_;
    typedef State<MODEL>                    State_;
    typedef State4D<MODEL>                  State4D_;

 public:
    CostPert(const eckit::Configuration &, const eckit::mpi::Comm &, const eckit::Configuration &,
             const bool & shiftTime = false);
    virtual ~CostPert();

    double evaluate(CtrlVar_ &, const eckit::Configuration &, PostProcessor<State_>) override;
    void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
                PostProcessor<Increment_>,
                const bool idModel = false) const override;
    void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
                PostProcessor<Increment_>,
                const bool idModel = false) const override;
    void zeroAD(CtrlInc_ &) const override;

    void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

    void resetLinearization() override {}

    const CostBase_ & jterm(const size_t ii) const override {return Jmem_->jterm(ii);}
    size_t nterms() const override {return Jmem_->nterms();}
    double getCostJoJc() const override {return Jmem_->getCostJoJc();}

 private:
    void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_> &) const override;

    int iteration_ = 0;
    const eckit::Configuration & contConf_;
    const util::TimeWindow timeWindow_;
    const eckit::Configuration & memConf_;
    const eckit::mpi::Comm & comm_;
    bool shiftTime_;
    std::unique_ptr<CostFct3DVar_> Jmem_;
    std::shared_ptr<CtrlVar_> xxMember_;
    State_ * hackBG_;
    std::shared_ptr<StateSaver<State_>> saver_;
};
// =============================================================================

template<typename MODEL, typename OBS, class J>
CostPert<MODEL, OBS, J>::CostPert(const eckit::Configuration & controlConf,
                                  const eckit::mpi::Comm & comm,
                                  const eckit::Configuration & memberConf,
                                  const bool & shiftTime)
  :J(controlConf.getSubConfiguration("cost function"), comm), contConf_(controlConf),
    timeWindow_(controlConf.getSubConfiguration("cost function.time window")),
    memConf_(memberConf), comm_(comm), shiftTime_(shiftTime), Jmem_(), xxMember_(),
    hackBG_(nullptr), saver_()
{
  Log::trace() << "CostPert::CostPert" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS, class J>
CostPert<MODEL, OBS, J>::~CostPert() {
  Jmem_.reset();
  Log::trace() << "CostPert::~CostPert" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS, class J>
double CostPert<MODEL, OBS, J>::evaluate(CtrlVar_ & fguess,
                                         const eckit::Configuration & innerConf,
                                         PostProcessor<State_> post) {
  if (iteration_ == 0) {
    PostProcessor<State_> postControl(post);
    std::vector<eckit::LocalConfiguration> iterconfs;
    contConf_.get("variational.iterations", iterconfs);
    iterconfs[0].set("linearize", true);
    iterconfs[0].set("iteration", 0);
    iterconfs[0].set("total iterations", innerConf.getInt("total iterations"));
    CtrlVar_ xxControl(this->jb().getBackground());
    if (shiftTime_) {
      hackBG_ = &xxControl.state();
      std::vector<std::string> antime = {timeWindow_.midpoint().toString()};
      eckit::LocalConfiguration halfwin;
      halfwin.set("times", antime);
      halfwin.set("variables", hackBG_->variables().variables());
      saver_.reset(new StateSaver<State_>(halfwin));
      postControl.enrollProcessor(saver_);
    }
    J::evaluate(fguess, iterconfs[0], postControl);

//  Set up the ensemble member background
//  Need shared pointers to State4D, ModelAuxControl & ObsAuxControl
    std::shared_ptr<ModelAuxControl_>
            modBias(new ModelAuxControl_(this->jb().jbState().geometry(),
                                         memConf_.getSubConfiguration("model aux error")));
    std::shared_ptr<ObsAuxControls_>
            obsBias(new ObsAuxControls_(this->jo().obspaces(),
                                        memConf_.getSubConfiguration("observations.observers")));
    std::shared_ptr<State4D_>
            state(new State4D_(this->jb().jbState().geometry(),
                               eckit::LocalConfiguration(memConf_, "background")));
    xxMember_.reset(new CtrlVar_(state, modBias, obsBias));

//  Change background to be the difference between the ensemble member and the control
    CtrlInc_ xxPertInc(this->jb());

    if (shiftTime_) {
      xxPertInc.state().updateTime(timeWindow_.length()/2);
      *hackBG_ = saver_->getState();    // the state component of xxControl is updated via hackBG_
    }

    xxPertInc.diff(*xxMember_, xxControl);
    CtrlVar_ xxPert(*xxMember_, false);
    this->addIncrement(xxPert, xxPertInc);

    CtrlVar_ xxPertCopy(xxPert);
    this->getNonConstJb()->getBackground() = xxPertCopy;
//  Set the obs operator to the linear obs operator
//  Perturb zero-valued obs ready for the pert member
    this->getNonConstJo()->setObsPert(xxPertInc.state().geometry(), xxPertInc.state().variables());

//  Set up pert member 3DVar cost function
//  Runs 3DVar methods using pointers to the Jb & Jo already set up
    Jmem_.reset(new CostFct3DVar_(memConf_, comm_, this->getNonConstJb(), this->getNonConstJo()));
    fguess = xxPert;
  }
//  Call evaluate on the pert cost function, without linearising
//  Iteration is set to more than 0 to ensure filters aren't rerun
  iteration_++;
  eckit::LocalConfiguration innerPertConf(innerConf);
  innerPertConf.set("iteration", iteration_);
  innerPertConf.set("linearize", false);
  innerPertConf.set("control pert", true);

  return Jmem_->evaluate(fguess, innerPertConf, post);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS, class J>
void CostPert<MODEL, OBS, J>::runTLM(CtrlInc_ & ctrlInc, PostProcessorTLAD<MODEL> & ppTLAD,
                                     PostProcessor<Increment_> post,
                                     const bool idModel) const {
  Log::trace() << "CostPert::runTLM start" << std::endl;
  Jmem_->runTLM(ctrlInc, ppTLAD, post, idModel);
  Log::trace() << "CostPert::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS, class J>
void CostPert<MODEL, OBS, J>::runADJ(CtrlInc_ & ctrlInc, PostProcessorTLAD<MODEL> & ppTLAD,
                                     PostProcessor<Increment_> post,
                                     const bool idModel) const {
  Log::trace() << "CostPert::runADJ start" << std::endl;
  Jmem_->runADJ(ctrlInc, ppTLAD, post, idModel);
  Log::trace() << "CostPert::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS, class J>
void CostPert<MODEL, OBS, J>::zeroAD(CtrlInc_ & ctrlInc) const {
  Log::trace() << "CostPert::zeroAD start" << std::endl;
  Jmem_->zeroAD(ctrlInc);
  Log::trace() << "CostPert::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS, class J>
void CostPert<MODEL, OBS, J>::runNL(CtrlVar_ & ctrlVar, PostProcessor<State_>& post) const {
  Log::trace() << "CostPert::runNL pert start" << std::endl;
  if (Jmem_ != nullptr && ctrlVar.state().validTime() != timeWindow_.start()) {
//  Jmem_ calls runNL for the pert assimilation or if the control cost function is 3DVar
    Jmem_->runNL(ctrlVar, post);
  } else {
//  J is used for the first call to evaluate,
//  and for the final analysis trajectory if the state is valid at the start of the window (4DVar)
    J::runNL(ctrlVar, post);
  }
  Log::trace() << "CostPert::runNL pert done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS, class J>
void CostPert<MODEL, OBS, J>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                      PostProcessor<Increment_>& post) const {
  Log::trace() << "CostPert::addIncr start" << std::endl;
  if (iteration_ == 0) {
    Log::trace() << "CostPert::addIncr doesn't exclude out of bounds values" << std::endl;
    atlas::FieldSet fset = xx.state().fieldSet().fieldSet();
    util::addFieldSets(fset, dx.state().fieldSet().fieldSet());
    xx.state().synchronizeFields();
  } else {
    Log::trace() << "CostPert::addIncr removes out of bounds values" << std::endl;
    ASSERT(dx.state().validTime() == timeWindow_.midpoint());
    // If the state passed in is valid at the start of the window,
    // but the increment is valid at the middle, the state is then set to be
    // equal to the mid window control member state. This works because the control
    // member is the only member that could have a state valid at the start of the window
    // and so avoids repeating the model integration. Note that the state will be valid at the
    // middle of the window following the call to addIncr.
    if (xx.state().validTime() == timeWindow_.start()) xx.state() = saver_->getState();
    ASSERT(xx.state().validTime() == timeWindow_.midpoint());
    xx.state() += dx.state();
  }
  Log::trace() << "CostPert::addIncr done" << std::endl;
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTPERT_H_
