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

template <typename MODEL>
class CalcHofX {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef State<MODEL>               State_;
  typedef State4D<MODEL>             State4D_;
  typedef PostProcessor<State_>      PostProcessor_;
  template <typename DATA> using ObsData_ = ObsDataVector<MODEL, DATA>;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsData_<DATA> >;

 public:
/// \brief Initializes Observers
  CalcHofX(const ObsSpaces_ &, const Geometry_ &, const eckit::LocalConfiguration &);

/// \brief Computes 4D H(x) (running the model)
  const Observations_ & compute(const Model_ &, State_ &, PostProcessor_ &);
/// \brief Computes 4D H(x) (using State4D)
  const Observations_ & compute(const State4D_ &);

/// \brief accessor to QC flag
  const ObsData_<int> & qcFlags(const size_t ii) const {return *(qcflags_[ii]);}
/// \brief accessor to Obs errors
  const ObsData_<float> & obsErrors(const size_t ii) const {return *(obserr_[ii]);}

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
  std::vector<ObsDataPtr_<int> > qcflags_;   // QC flags
  std::vector<ObsDataPtr_<float> > obserr_;  // Obs Errors
  boost::shared_ptr<Observers<MODEL, State_> > pobs_;   // Observer
};


// -----------------------------------------------------------------------------

template <typename MODEL>
CalcHofX<MODEL>::CalcHofX(const ObsSpaces_ & obspaces, const Geometry_ & geometry,
                          const eckit::LocalConfiguration & config) :
  obsconf_(config.getSubConfiguration("Observations")),
  obspaces_(obspaces), geometry_(geometry),
  ybias_(obspaces_, config.getSubConfiguration("Observations")),
  moderr_(geometry_, config.getSubConfiguration("Initial Condition")),
  winbgn_(config.getString("Assimilation Window.window_begin")),
  winlen_(config.getString("Assimilation Window.window_length")) {}

// -----------------------------------------------------------------------------
template <typename MODEL>
void CalcHofX<MODEL>::initObserver() {
  qcflags_.clear();
  obserr_.clear();
  qcflags_.reserve(obspaces_.size());
  obserr_.reserve(obspaces_.size());
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
//  Allocate QC flags
    qcflags_.emplace_back(boost::make_shared<ObsData_<int>>(obspaces_[jj],
                            obspaces_[jj].obsvariables()));
//  Allocate and read initial obs error
    obserr_.emplace_back(boost::make_shared<ObsData_<float>>(obspaces_[jj],
                            obspaces_[jj].obsvariables(), "ObsError"));
  }
//  Setup Observers
  pobs_.reset(new Observers<MODEL, State_>(obsconf_, obspaces_, ybias_, qcflags_, obserr_));
  oops::Log::trace() << "CalcHofX<MODEL>::CalcHofX created" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
const Observations<MODEL> & CalcHofX<MODEL>::compute(const Model_ & model, State_ & xx,
                                                     PostProcessor_ & post) {
  oops::Log::trace() << "CalcHofX<MODEL>::compute (model) start" << std::endl;

  this->initObserver();
//  run the model and compute H(x)
  post.enrollProcessor(pobs_);
  model.forecast(xx, moderr_, winlen_, post);

  oops::Log::trace() << "CalcHofX<MODEL>::compute (model) done" << std::endl;
  return pobs_->hofx();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
const Observations<MODEL> & CalcHofX<MODEL>::compute(const State4D_ & xx) {
  oops::Log::trace() << "CalcHofX<MODEL>::compute (state4D) start" << std::endl;

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

  oops::Log::trace() << "CalcHofX<MODEL>::compute (state4D) done" << std::endl;
  return pobs_->hofx();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_CALCHOFX_H_
