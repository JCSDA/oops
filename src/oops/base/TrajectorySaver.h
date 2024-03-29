/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_TRAJECTORYSAVER_H_
#define OOPS_BASE_TRAJECTORYSAVER_H_

#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/LinearModel.h"
#include "oops/base/PostBase.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/PostTimerParameters.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class TrajectorySaver4DParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(TrajectorySaver4DParameters, Parameters)

 public:
  RequiredPolymorphicParameter<LinearModelParametersBase, LinearModelFactory<MODEL>>
    tlmParameters{"name", this};
  PostTimerParameters postTimer{this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class TrajectorySaver3DParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(TrajectorySaver3DParameters, Parameters)

 public:
  PostTimerParameters postTimer{this};
  // Some places in the code still pass Configuration objects to constructors of TrajectorySaver.
  // These objects often contain settings unrelated to TrajectorySaver, which would normally be
  // flagged as errors (potential misspellings) when the contents of these objects are validated
  // against the JSON schema defined by TrajectorySaver3DParameters. The ignoreOthers member causes
  // the schema to ignore unknown top-level settings.
  //
  // This problem will disappear when these parts of the code are refactored to define appropriate
  // Parameters subclasses, typically containing a member of type TrajectorySaver3DParameters that
  // will be passed to the constructor of TrajectorySaver instead of the Configuration object.
  //
  // TODO(wsmigaj): Remove ignoreOthers (and the TrajectorySaver constructors taking Configuration
  // objects) when it becomes possible.
  IgnoreOtherParameters ignoreOthers{this};
};

/// Save trajectory during forecast run.

// -----------------------------------------------------------------------------

template <typename MODEL>
class TrajectorySaver : public PostBase<State<MODEL> > {
  typedef Geometry<MODEL>          Geometry_;
  typedef LinearModel<MODEL>       LinearModel_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef PostProcessorTLAD<MODEL> PPTLAD_;
  typedef State<MODEL>             State_;
  typedef TrajectorySaver4DParameters<MODEL> Parameters4D_;
  typedef TrajectorySaver3DParameters<MODEL> Parameters3D_;

 public:
  TrajectorySaver(const Parameters4D_ &, const Geometry_ &,
                  const ModelAux_ &, std::shared_ptr<LinearModel_>, PPTLAD_);
  TrajectorySaver(const eckit::Configuration &, const Geometry_ &,
                  const ModelAux_ &, std::shared_ptr<LinearModel_>, PPTLAD_);
  TrajectorySaver(const Parameters3D_ &, const Geometry_ &, PPTLAD_);
  TrajectorySaver(const eckit::Configuration &, const Geometry_ &, PPTLAD_);
  ~TrajectorySaver() {}

 private:
  const Geometry_ & resol_;
  PPTLAD_ pptraj_;
  std::unique_ptr<const ModelAux_>  lrBias_;
  std::shared_ptr<LinearModel_>     tlm_;

  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
TrajectorySaver<MODEL>::TrajectorySaver(const Parameters4D_ & parameters,
                                        const Geometry_ & resol,
                                        const ModelAux_ & bias,
                                        std::shared_ptr<LinearModel_> tlm,
                                        PPTLAD_ pptraj):
  PostBase<State_>(parameters.postTimer),
  resol_(resol), pptraj_(pptraj), lrBias_(new ModelAux_(resol, bias)), tlm_(tlm)
{
  Log::trace() << "TrajectorySaver::TrajectorySaver 4D" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
TrajectorySaver<MODEL>::TrajectorySaver(const eckit::Configuration & conf,
                                        const Geometry_ & resol,
                                        const ModelAux_ & bias,
                                        std::shared_ptr<LinearModel_> tlm,
                                        PPTLAD_ pptraj):
  TrajectorySaver(validateAndDeserialize<Parameters4D_>(conf), resol, bias, tlm, pptraj)
{}
// -----------------------------------------------------------------------------
template <typename MODEL>
TrajectorySaver<MODEL>::TrajectorySaver(const Parameters3D_ & parameters,
                                        const Geometry_ & resol, PPTLAD_ pptraj):
  PostBase<State_>(parameters.postTimer),
  resol_(resol), pptraj_(pptraj), lrBias_(), tlm_()
{
  Log::trace() << "TrajectorySaver::TrajectorySaver 3D" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
TrajectorySaver<MODEL>::TrajectorySaver(const eckit::Configuration & conf,
                                        const Geometry_ & resol, PPTLAD_ pptraj):
  TrajectorySaver(validateAndDeserialize<Parameters3D_>(conf), resol, pptraj)
{}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doInitialize(const State_ & x0,
                                          const util::DateTime & end,
                                          const util::Duration & step) {
  Log::trace() << "TrajectorySaver::doInitialize start" << std::endl;
  State_ xlr(resol_, x0);
  pptraj_.initializeTraj(xlr, end, step);
  Log::trace() << "TrajectorySaver::doInitialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doProcessing(const State_ & xx) {
  Log::trace() << "TrajectorySaver::doProcessing start" << std::endl;
  State_ xlr(resol_, xx);
  if (tlm_) tlm_->setTrajectory(xx, xlr, *lrBias_);
  pptraj_.processTraj(xlr);
  Log::trace() << "TrajectorySaver::doProcessing done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doFinalize(const State_ & xx) {
  Log::trace() << "TrajectorySaver::doFinalize start" << std::endl;
  State_ xlr(resol_, xx);
  pptraj_.finalizeTraj(xlr);
  Log::trace() << "TrajectorySaver::doFinalize done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_TRAJECTORYSAVER_H_
