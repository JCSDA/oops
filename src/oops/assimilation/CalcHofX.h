/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_ASSIMILATION_CALCHOFX_H_
#define OOPS_ASSIMILATION_CALCHOFX_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "oops/assimilation/State4D.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/QCData.h"
#include "oops/base/StateInfo.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Computes observation operator (while running model, or with State4D)

template <typename MODEL, typename OBS>
class CalcHofX {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<OBS>        ObsAuxCtrls_;
  typedef Observations<OBS>          Observations_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State<MODEL>               State_;
  typedef State4D<MODEL>             State4D_;
  typedef PostProcessor<State_>      PostProcessor_;
  typedef QCData<OBS>                QCData_;
  template <typename DATA> using ObsData_ = ObsDataVector<OBS, DATA>;

 public:
/// \brief Initializes Observers
  CalcHofX(const ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &);

/// \brief Computes 4D H(x) (running the model)
  const Observations_ & compute(const Model_ &, State_ &, PostProcessor_ &);
/// \brief Computes 4D H(x) (using State4D)
  const Observations_ & compute(const State4D_ &);

/// \brief saves QC flags to ObsSpaces
  void saveQcFlags(const std::string &) const;
/// \brief saves obs error variances (modified in QC) to ObsSpaces
  void saveObsErrors(const std::string &) const;

 private:
/// \brief helper method to initialize qc flags and observer
  void initObserver();

  const eckit::LocalConfiguration obsconf_;  // configuration for observer
  const ObsSpaces_ & obspaces_;              // ObsSpaces used in H(x)
  const Geometry_ &  geometry_;              // Model Geometry
  ObsAuxCtrls_       ybias_;                 // obs bias
  ModelAux_          moderr_;                // model bias
  const util::DateTime winbgn_;              // window for assimilation
  const util::Duration winlen_;
  std::unique_ptr<QCData_> qc_;              // QC-related (flags and obserrors)
  std::shared_ptr<Observers<MODEL, OBS> > pobs_;  // Observer
};


// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CalcHofX<MODEL, OBS>::CalcHofX(const ObsSpaces_ & obspaces, const Geometry_ & geometry,
                          const eckit::Configuration & config) :
  obsconf_(config.getSubConfiguration("Observations")),
  obspaces_(obspaces), geometry_(geometry),
  ybias_(obspaces_, config.getSubConfiguration("Observations")),
  moderr_(geometry_, config.getSubConfiguration("Initial Condition")),
  winbgn_(config.getString("Assimilation Window.window_begin")),
  winlen_(config.getString("Assimilation Window.window_length")) {}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void CalcHofX<MODEL, OBS>::initObserver() {
  qc_.reset(new QCData_(obspaces_));
//  Setup Observers
  pobs_.reset(new Observers<MODEL, OBS>(obsconf_, obspaces_, ybias_, *qc_));
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
const Observations<OBS> & CalcHofX<MODEL, OBS>::compute(const Model_ & model, State_ & xx,
                                                     PostProcessor_ & post) {
  oops::Log::trace() << "CalcHofX<MODEL, OBS>::compute (model) start" << std::endl;

  this->initObserver();
//  run the model and compute H(x)
  post.enrollProcessor(pobs_);
  model.forecast(xx, moderr_, winlen_, post);

  oops::Log::trace() << "CalcHofX<MODEL, OBS>::compute (model) done" << std::endl;
  return pobs_->hofx();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
const Observations<OBS> & CalcHofX<MODEL, OBS>::compute(const State4D_ & xx) {
  oops::Log::trace() << "CalcHofX<MODEL, OBS>::compute (state4D) start" << std::endl;

  this->initObserver();
  size_t nstates = xx.size();
  util::DateTime winend = winbgn_ + winlen_;
  util::Duration tstep = winlen_;   // for a single state
  // if using several states, compute the timestep and check that it's the same
  // for all states
  if (nstates > 1) {
    tstep = xx[1].validTime() - xx[0].validTime();
    for (size_t ii = 1; ii < xx.size(); ++ii) {
      ASSERT(tstep == (xx[ii].validTime() - xx[ii-1].validTime()));
    }
  }
  // check that initial and last states have valid times
  ASSERT(xx[0].validTime() <= (winbgn_ + tstep/2));
  ASSERT(xx[nstates-1].validTime() >= (winend - tstep/2));

  // run Observer looping through all the states
  pobs_->initialize(xx[0], winend, tstep);
  for (size_t ii = 0; ii < xx.size(); ++ii) {
    pobs_->process(xx[ii]);
  }
  pobs_->finalize(xx[nstates-1]);

  oops::Log::trace() << "CalcHofX<MODEL, OBS>::compute (state4D) done" << std::endl;
  return pobs_->hofx();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CalcHofX<MODEL, OBS>::saveQcFlags(const std::string & name) const {
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    qc_->qcFlags(jj)->save(name);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CalcHofX<MODEL, OBS>::saveObsErrors(const std::string & name) const {
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    qc_->obsErrors(jj)->save(name);
  }
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_CALCHOFX_H_
