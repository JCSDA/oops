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
#include <numeric>
#include <string>
#include <vector>

#include "oops/base/Increment4D.h"
#include "oops/base/Model.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/generic/SimpleLinearModel.h"

namespace oops {

/*
 * Configuration options for HtlmEnsemble:
 *
 * Top-level keys:
 * ─────────────────────────────────────────────────────────────────────────────
 * "model"                  : Configuration for the model used in nonlinear control.
 * "model geometry"         : Configuration for the geometry used in the nonlinear control.
 * "model for ensemble"     : (Optional) Configuration for the model used in ensemble members.
 * "geometry for ensemble"  : (Optional) Geometry for ensemble; defaults to control geometry.
 * "nonlinear control"      : Configuration for the nonlinear control/trajectory initial condition.
 * "nonlinear ensemble"     : Configuration block for initializing the nonlinear ensemble.
 *
 * Subkeys under "nonlinear ensemble":
 * ─────────────────────────────────────────────────────────────────────────────
 * "read"                   : (Optional) Read ensemble initial conditions from file.
 *
 * "generate"               : (Optional) Generate ensemble using covariance perturbations.
 *     → contains:
 *        - "ensemble size"      : Number of ensemble members to generate.
 *        - "background error"   : Configuration for the background error covariance.
 *        - "variables"          : List of variables to perturb.
 *
 * Notes:
 * - Either "read" or "generate" must be provided under "nonlinear ensemble", but not both.
 * - If both are present, the constructor will abort with an error.
 */

//------------------------------------------------------------------------------

template <typename MODEL>
class HtlmEnsemble{
  typedef CovarianceFactory<MODEL>                     CovarianceFactory_;
  typedef Geometry<MODEL>                              Geometry_;
  typedef Increment<MODEL>                             Increment_;
  typedef Increment4D<MODEL>                           Increment4D_;
  typedef IncrementSet<MODEL>                          IncrementSet_;
  typedef Model<MODEL>                                 Model_;
  typedef ModelAuxControl<MODEL>                       ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>                     ModelAuxIncrement_;
  typedef ModelSpaceCovarianceBase<MODEL>              CovarianceBase_;
  typedef SimpleLinearModel<MODEL>                     SimpleLinearModel_;
  typedef State<MODEL>                                 State_;
  typedef StateSet<MODEL>                              StateSet_;
  typedef State4D<MODEL>                               State4D_;

 public:
  static const std::string classname() {return "oops::HtlmEnsemble";}

  HtlmEnsemble(const eckit::Configuration &, SimpleLinearModel_ &, const Geometry_ &,
               const Variables &);
  atlas::FieldSet getRmsVals(const Variables &, const atlas::idx_t) const;
  void step(const util::Duration &, SimpleLinearModel_ &);

  IncrementSet_ & getLinearEnsemble() {return linearEnsemble_;}
  const IncrementSet_ & getLinearEnsemble() const {return linearEnsemble_;}
  const IncrementSet_ & getLinearErrors() const {return linearErrors_;}
  size_t size() const {return ensembleSize_;}

 private:
  static const std::vector<int> getEnsVec(const size_t & ensSize) {
    std::vector<int> ensVec(ensSize);
    std::iota(ensVec.begin(), ensVec.end(), 0);
    return ensVec;
  }
  const eckit::LocalConfiguration nlEnsConf_;
  const eckit::LocalConfiguration ensModelConf_;

  const Geometry_ & updateGeometry_;
  std::shared_ptr<Geometry_> controlGeometry_;
  std::shared_ptr<Geometry_> ensembleGeometry_;
  std::shared_ptr<Model_> modelControl_;
  std::shared_ptr<Model_> modelEnsemble_;
  State4D_ nonlinearControl_;
  StateSet_ nonlinearEnsemble_;
  std::unique_ptr<State_> spareStateEnsembleGeometry_;
  const size_t ensembleSize_;
  IncrementSet_ nonlinearDifferences_;
  IncrementSet_ linearEnsemble_;
  IncrementSet_ linearErrors_;
  ModelAuxCtl_ maux_;
  ModelAuxIncrement_ mauxinc_;
  PostProcessor<State_> trajectorySaver_;
  PostProcessor<State_> emptyPp_;
};

//------------------------------------------------------------------------------

template<typename MODEL>
HtlmEnsemble<MODEL>::HtlmEnsemble(const eckit::Configuration & config,
                                  SimpleLinearModel_ & simpleLinearModel,
                                  const Geometry_ & updateGeometry,
                                  const Variables & vars)
: nlEnsConf_(config, "nonlinear ensemble"), updateGeometry_(updateGeometry),
  controlGeometry_(std::make_shared<Geometry_>(
    config.getSubConfiguration("model geometry"), updateGeometry_.getComm())),
  ensembleGeometry_(config.has("geometry for ensemble")
    ? std::make_shared<Geometry_>(config.getSubConfiguration("geometry for ensemble"),
      updateGeometry_.getComm())
    : controlGeometry_),
  modelControl_(std::make_shared<Model_>(
    *controlGeometry_, config.getSubConfiguration("model"))),
  modelEnsemble_(config.has("model for ensemble") ? std::make_shared<Model_>(
    *ensembleGeometry_, config.getSubConfiguration("model for ensemble"))
    : modelControl_),
  nonlinearControl_(*controlGeometry_, config.getSubConfiguration("nonlinear control")),
  nonlinearEnsemble_(nlEnsConf_.has("read") ?
    StateSet_(*ensembleGeometry_,
      nlEnsConf_.getSubConfiguration("read")) :
    StateSet_(*ensembleGeometry_, nonlinearControl_[0].variables(),
      {nonlinearControl_[0].validTime()}, oops::mpi::myself(),
        getEnsVec(nlEnsConf_.getSubConfiguration("generate")
          .getUnsigned("ensemble size")))),
  spareStateEnsembleGeometry_(controlGeometry_ == ensembleGeometry_ ?
    nullptr : std::make_unique<State_>(*ensembleGeometry_, nonlinearControl_[0])),
  ensembleSize_(nonlinearEnsemble_.size()),
  nonlinearDifferences_(*ensembleGeometry_, vars, {nonlinearControl_[0].validTime()},
                        oops::mpi::myself(), getEnsVec(ensembleSize_)),
  linearEnsemble_(updateGeometry_, vars, {nonlinearControl_[0].validTime()},
                        oops::mpi::myself(), getEnsVec(ensembleSize_)),
  linearErrors_(linearEnsemble_), maux_(*ensembleGeometry_, eckit::LocalConfiguration()),
  mauxinc_(updateGeometry_, eckit::LocalConfiguration())
{
  Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() starting" << std::endl;
  // If required, initialize nonlinear ensemble from covariance
  if (nlEnsConf_.has("generate")) {
    if (nlEnsConf_.has("read")) {
      ABORT("HtlmEnsemble<MODEL>: both types of nonlinear ensemble initial conditions provided");
    }
    const eckit::LocalConfiguration genConf(config, "nonlinear ensemble.generate");
    const Variables genVars(genConf, "variables");
    const eckit::LocalConfiguration covConf(genConf, "background error");
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
      *ensembleGeometry_, genVars, covConf, nonlinearControl_, nonlinearControl_));
    Increment4D_ dx(*ensembleGeometry_, genVars, nonlinearControl_.times());
    for (size_t m = 0; m < ensembleSize_; m++) {
      Bmat->randomize(dx);
      nonlinearEnsemble_[m] = nonlinearControl_[0];
      nonlinearEnsemble_[m] += dx[0];
    }
  }
  // Set up linearEnsemble_ initial conditions
  Increment_ linearEnsembleMemberEnsembleGeometry(*ensembleGeometry_, vars,
                                                  nonlinearControl_[0].validTime());
  for (size_t m = 0; m < ensembleSize_; m++) {
    linearEnsembleMemberEnsembleGeometry.diff(nonlinearEnsemble_[m],
      controlGeometry_ == ensembleGeometry_ ? nonlinearControl_[0] : *spareStateEnsembleGeometry_);
    linearEnsemble_[m] = Increment_(updateGeometry_, linearEnsembleMemberEnsembleGeometry);
  }
  // Set up a TrajectorySaver for simpleLinearModel_
  simpleLinearModel.setUpTrajectorySaver(trajectorySaver_, maux_);
  Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
atlas::FieldSet HtlmEnsemble<MODEL>::getRmsVals(const Variables & updateVars,
                                                const atlas::idx_t nLevels) const {
  atlas::FieldSet rmsVals;
  for (const auto & var : updateVars) {
    std::vector<double> rmsVar = linearEnsemble_[0].rmsByVariableByLevel(var, false);
    rmsVals.add(atlas::Field(
      var.name(), atlas::array::make_datatype<double>(), atlas::array::make_shape(nLevels)));
    auto rmsView = atlas::array::make_view<double, 1>(rmsVals[var.name()]);
    // Avoid divide-by-zero
    for (atlas::idx_t k = 0; k < nLevels; k++) rmsView[k] = (rmsVar[k]) ? rmsVar[k] : 1.0;
  }
  return rmsVals;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmEnsemble<MODEL>::step(const util::Duration & tstep,
                               SimpleLinearModel_ & simpleLinearModel) {
  Log::trace() << "HtlmEnsemble<MODEL>::step() starting" << std::endl;
  modelControl_->forecast(nonlinearControl_[0], maux_, tstep, trajectorySaver_);
  if (controlGeometry_ != ensembleGeometry_) {
    *spareStateEnsembleGeometry_ = State_(*ensembleGeometry_, nonlinearControl_[0]);
  }
  for (size_t m = 0; m < ensembleSize_; m++) {
    modelEnsemble_->forecast(nonlinearEnsemble_[m], maux_, tstep, emptyPp_);
    nonlinearDifferences_[m].updateTime(tstep);
    nonlinearDifferences_[m].diff(nonlinearEnsemble_[m],
      controlGeometry_ == ensembleGeometry_ ? nonlinearControl_[0] : *spareStateEnsembleGeometry_);
    linearErrors_[m] = Increment_(updateGeometry_, nonlinearDifferences_[m]);
    simpleLinearModel.forecastTL(linearEnsemble_[m], mauxinc_, tstep);
    linearErrors_[m] -= linearEnsemble_[m];
  }
  Log::trace() << "HtlmEnsemble<MODEL>::step() done" << std::endl;
}

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMENSEMBLE_H_
