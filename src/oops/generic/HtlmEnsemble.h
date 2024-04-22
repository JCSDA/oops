/*
 * (C) Copyright 2022-2023 UCAR.
 * (C) Crown copyright 2022-2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMENSEMBLE_H_
#define OOPS_GENERIC_HTLMENSEMBLE_H_

#include <memory>
#include <string>

#include "oops/base/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Model.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/generic/SimpleLinearModel.h"

namespace oops {

template <typename MODEL>
class StatePerturbationParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StatePerturbationParameters, Parameters)
  typedef ModelSpaceCovarianceParametersWrapper<MODEL>    CovarianceParameters_;

 public:
  RequiredParameter<size_t> ensembleSize{"ensemble size", this};
  RequiredParameter<CovarianceParameters_> backgroundError{"background error", this};
  RequiredParameter<Variables> variables{"variables", this};
};

//------------------------------------------------------------------------------

template<typename MODEL>
class NonlinearEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(NonlinearEnsembleParameters, Parameters)
  typedef StatePerturbationParameters<MODEL>    StatePerturbationParameters_;
  typedef StateEnsembleParameters<MODEL>        StateEnsembleParameters_;

 public:
  OptionalParameter<StateEnsembleParameters_> fromFile{
    "read", "read nonlinear ensemble initial conditions from file", this};
  OptionalParameter<StatePerturbationParameters_> fromCovar{
    "generate", "generate nonlinear ensemble initial conditions from covariance", this};

  void check() const {
    if (fromFile.value() == boost::none && fromCovar.value() == boost::none) {
      ABORT("NonlinearEnsembleParameters<MODEL>: "
            "no nonlinear ensemble initial conditions provided (HtlmEnsemble)");
    }
    if (fromFile.value() != boost::none && fromCovar.value() != boost::none) {
      ABORT("NonlinearEnsembleParameters<MODEL>: "
            "both types of nonlinear ensemble initial conditions provided (HtlmEnsemble)");
    }
  }
};

//------------------------------------------------------------------------------

template <typename MODEL>
class HtlmEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmEnsembleParameters, Parameters);
  typedef typename Geometry<MODEL>::Parameters_    GeometryParameters_;
  typedef NonlinearEnsembleParameters<MODEL>       NonlinearEnsembleParameters_;
  typedef StateEnsembleParameters<MODEL>           StateEnsembleParameters_;

 public:
  RequiredParameter<eckit::LocalConfiguration> model{"model", this};
  RequiredParameter<GeometryParameters_> modelGeometry{"model geometry", this};
  RequiredParameter<eckit::LocalConfiguration> nonlinearControl{"nonlinear control", this};
  RequiredParameter<NonlinearEnsembleParameters_> nonlinearEnsemble{"nonlinear ensemble", this};
};

//------------------------------------------------------------------------------

template <typename MODEL>
class HtlmEnsemble{
  typedef CovarianceFactory<MODEL>                     CovarianceFactory_;
  typedef Geometry<MODEL>                              Geometry_;
  typedef HtlmEnsembleParameters<MODEL>                Parameters_;
  typedef Increment<MODEL>                             Increment_;
  typedef Increment4D<MODEL>                           Increment4D_;
  typedef IncrementEnsemble<MODEL>                     IncrementEnsemble_;
  typedef Model<MODEL>                                 Model_;
  typedef ModelAuxControl<MODEL>                       ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>                     ModelAuxIncrement_;
  typedef ModelSpaceCovarianceBase<MODEL>              CovarianceBase_;
  typedef ModelSpaceCovarianceParametersBase<MODEL>    CovarianceParameters_;
  typedef SimpleLinearModel<MODEL>                     SimpleLinearModel_;
  typedef State<MODEL>                                 State_;
  typedef StateEnsemble<MODEL>                         StateEnsemble_;
  typedef StatePerturbationParameters<MODEL>           PerturbationParameters_;
  typedef State4D<MODEL>                               State4D_;

 public:
  static const std::string classname() {return "oops::HtlmEnsemble";}

  HtlmEnsemble(const Parameters_ &, SimpleLinearModel_ &, const Geometry_ &, const Variables &);
  void step(const util::Duration &, SimpleLinearModel_ &);

  IncrementEnsemble_ & getLinearEnsemble() {return linearEnsemble_;}
  const IncrementEnsemble_ & getLinearEnsemble() const {return linearEnsemble_;}
  const IncrementEnsemble_ & getLinearErrors() const {return linearErrors_;}
  size_t size() const {return ensembleSize_;}

 private:
  const Geometry_ & updateGeometry_;
  const Geometry_ modelGeometry_;
  const Model_ model_;
  State4D_ nonlinearControl_;
  StateEnsemble_ nonlinearEnsemble_;
  const size_t ensembleSize_;
  IncrementEnsemble_ nonlinearDifferences_;
  IncrementEnsemble_ linearEnsemble_;
  IncrementEnsemble_ linearErrors_;
  ModelAuxCtl_ maux_;
  ModelAuxIncrement_ mauxinc_;
  PostProcessor<State_> trajectorySaver_;
  PostProcessor<State_> emptyPp_;
};

//------------------------------------------------------------------------------

template<typename MODEL>
HtlmEnsemble<MODEL>::HtlmEnsemble(const Parameters_ & params,
                                  SimpleLinearModel_ & simpleLinearModel,
                                  const Geometry_ & updateGeometry,
                                  const Variables & vars)
: updateGeometry_(updateGeometry),
  modelGeometry_(params.modelGeometry.value(), updateGeometry_.getComm()),
  model_(modelGeometry_, eckit::LocalConfiguration(params.toConfiguration(), "model")),
  nonlinearControl_(modelGeometry_, params.nonlinearControl.value()),
  nonlinearEnsemble_(params.nonlinearEnsemble.value().fromFile.value() != boost::none ?
    StateEnsemble_(modelGeometry_, *params.nonlinearEnsemble.value().fromFile.value()) :
    StateEnsemble_(nonlinearControl_[0],
                   (*params.nonlinearEnsemble.value().fromCovar.value()).ensembleSize.value())),
  ensembleSize_(nonlinearEnsemble_.size()),
  nonlinearDifferences_(modelGeometry_, vars,  nonlinearControl_[0].validTime(), ensembleSize_),
  linearEnsemble_(updateGeometry_, vars, nonlinearControl_[0].validTime(), ensembleSize_),
  linearErrors_(linearEnsemble_), maux_(modelGeometry_, eckit::LocalConfiguration()),
  mauxinc_(updateGeometry_, eckit::LocalConfiguration())
{
  Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() starting" << std::endl;
  // If required, initialize nonlinear ensemble from covariance
  params.nonlinearEnsemble.value().check();
  if (params.nonlinearEnsemble.value().fromCovar.value() != boost::none) {
    const PerturbationParameters_ pertParams = *params.nonlinearEnsemble.value().fromCovar.value();
    const Variables vars(pertParams.variables);
    const CovarianceParameters_ & covParams
      = pertParams.backgroundError.value().covarianceParameters;
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
      modelGeometry_, vars, covParams, nonlinearControl_, nonlinearControl_));
    Increment4D_ dx(modelGeometry_, vars, nonlinearControl_.times());
    for (size_t m = 0; m < ensembleSize_; m++) {
      Bmat->randomize(dx);
      nonlinearEnsemble_[m] += dx[0];
    }
  }
  // Set up linearEnsemble_ initial conditions
  Increment_ linearEnsembleMemberModelGeometry(modelGeometry_, vars,
                                               nonlinearControl_[0].validTime());
  for (size_t m = 0; m < ensembleSize_; m++) {
    linearEnsembleMemberModelGeometry.diff(nonlinearControl_[0], nonlinearEnsemble_[m]);
    linearEnsemble_[m] = Increment_(updateGeometry_, linearEnsembleMemberModelGeometry);
  }
  // Set up a TrajectorySaver for simpleLinearModel_
  simpleLinearModel.setUpTrajectorySaver(trajectorySaver_, maux_);
  Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmEnsemble<MODEL>::step(const util::Duration & tstep,
                               SimpleLinearModel_ & simpleLinearModel) {
  Log::trace() << "HtlmEnsemble<MODEL>::step() starting" << std::endl;
  model_.forecast(nonlinearControl_[0], maux_, tstep, trajectorySaver_);
  for (size_t m = 0; m < ensembleSize_; m++) {
    model_.forecast(nonlinearEnsemble_[m], maux_, tstep, emptyPp_);
    nonlinearDifferences_[m].updateTime(tstep);
    nonlinearDifferences_[m].diff(nonlinearControl_[0], nonlinearEnsemble_[m]);
    linearErrors_[m] = Increment_(updateGeometry_, nonlinearDifferences_[m]);
    simpleLinearModel.forecastTL(linearEnsemble_[m], mauxinc_, tstep);
    linearErrors_[m] -= linearEnsemble_[m];
  }
  Log::trace() << "HtlmEnsemble<MODEL>::step() done" << std::endl;
}

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMENSEMBLE_H_
