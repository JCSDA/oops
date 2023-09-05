/*
 * (C) Crown Copyright 2023, the Met Office.
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

 private:
    void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_> &) const override;

    int iteration_ = 0;
    const eckit::Configuration & contConf_;
    const eckit::Configuration & memConf_;
    const eckit::mpi::Comm & comm_;
    bool shiftTime_;
    util::Duration windowLength_;
    util::DateTime windowBegin_;
    std::unique_ptr<CostFct3DVar_> Jmem_;
    std::shared_ptr<CtrlVar_> xxMember_;
};
// =============================================================================

template<typename MODEL, typename OBS, class J>
CostPert<MODEL, OBS, J>::CostPert(const eckit::Configuration & controlConf,
                                  const eckit::mpi::Comm & comm,
                                  const eckit::Configuration & memberConf,
                                  const bool & shiftTime)
  :J(controlConf.getSubConfiguration("cost function"), comm), contConf_(controlConf),
    memConf_(memberConf), comm_(comm), shiftTime_(shiftTime), Jmem_(), xxMember_()
{
  windowLength_ = util::Duration(
              controlConf.getSubConfiguration("cost function").getString("window length"));
  windowBegin_ = util::DateTime(
              controlConf.getSubConfiguration("cost function").getString("window begin"));
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
    std::vector<eckit::LocalConfiguration> iterconfs;
    contConf_.get("variational.iterations", iterconfs);
    iterconfs[0].set("linearize", true);
    iterconfs[0].set("iteration", 0);
    J::evaluate(fguess, iterconfs[0], post);

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
    xxPertInc.diff(*xxMember_, this->jb().getBackground());
    CtrlVar_ xxPert(*xxMember_, false);
    this->addIncrement(xxPert, xxPertInc);

    if (shiftTime_) {
        xxPert.state().updateTime(windowLength_/2);
    }

    CtrlVar_ xxPertCopy(xxPert);
    this->getNonConstJb()->getBackground() = xxPertCopy;
//  Set the obs operator to the linear obs operator
//  Perturb zero-valued obs ready for the pert member
    this->getNonConstJo()->setObsPert();

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
  if (Jmem_ != nullptr && ctrlVar.state().validTime() != windowBegin_) {
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
  if (iteration_ <= 1) {
    Log::trace() << "CostPert::addIncr doesn't exclude out of bounds values" << std::endl;
    atlas::FieldSet fset = xx.state().fieldSet();
    util::addFieldSets(fset, dx.state().fieldSet());
    xx.state().fromFieldSet(fset);
    xx.state().synchronizeFields();
  } else {
    Log::trace() << "CostPert::addIncr removes out of bounds values" << std::endl;
    xx.state() += dx.state();
  }
  Log::trace() << "CostPert::addIncr done" << std::endl;
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTPERT_H_
