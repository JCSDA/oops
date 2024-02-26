/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/PostBase.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class StateSetSaver : public PostBase<State<MODEL> > {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef StateSet<MODEL>            StateSet_;

 public:
  StateSetSaver(const eckit::Configuration &, const Geometry_ &);
  StateSetSaver(const eckit::Configuration &, const Geometry_ &,
                           const std::vector<util::DateTime> &,
                           const eckit::mpi::Comm &,
                           const std::vector<int> &,
                           const eckit::mpi::Comm &);
  ~StateSetSaver() {}
  void addState(const State_ & xx);
  std::unique_ptr<StateSet<MODEL> > & getStateSet(void);

 private:
  const Geometry_ & resol_;
  const std::vector<util::DateTime> & times_;
  const eckit::mpi::Comm & commTime_;
  const std::vector<int> & ens_;
  const eckit::mpi::Comm & commEns_;
  std::unique_ptr<StateSet_ > States_;
  size_t stateIndex_ = 0;
  bool initialized_ = false;

//  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doInitialize(const State_ &,
                    const util::DateTime &,
                    const util::Duration &) override;

  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
StateSetSaver<MODEL>::StateSetSaver(const eckit::Configuration & conf, const Geometry_ & resol):
  PostBase<State_>(conf), resol_(resol)
{
  Log::info() << "StateSetSaver::constructor" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
StateSetSaver<MODEL>::StateSetSaver(const eckit::Configuration & conf,
                           const Geometry_ & resol,
                           const std::vector<util::DateTime> & times,
                           const eckit::mpi::Comm & commTime,
                           const std::vector<int> & ens,
                           const eckit::mpi::Comm & commEns):
  PostBase<State_>(conf),
  resol_(resol),
  commTime_(commTime),
  commEns_(commEns),
  times_(times),
  ens_(ens)
{
  Log::info() << "StateSetSaver::constructor" << std::endl;
}
// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<StateSet<MODEL> > & StateSetSaver<MODEL>::getStateSet(void) {
  Log::info() << "StateSetSaver::returning StateSet" << std::endl;
  return(States_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void StateSetSaver<MODEL>::doInitialize(const State_ & x0,
                           const util::DateTime & bgndate,
                           const util::Duration & fcstlen ) {
  if (!initialized_) {
    Log::info() << "StateSetSaver::doInitialize start" << std::endl;
    States_.reset( new StateSet(resol_, x0.variables(), times_, commTime_, ens_, commEns_) );
    Log::info() << "StateSetSaver::doInitialize done" << std::endl;
    initialized_ = true;
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void StateSetSaver<MODEL>::doProcessing(const State_ & xx) {
  Log::info() << "StateSetSaver::doProcessing on stateIndex_ " << stateIndex_ << std::endl;
  (*States_)[stateIndex_] = xx;
  stateIndex_++;
  Log::info() << "StateSetSaver::doProcessing done" << std::endl;
}
template <typename MODEL>
void StateSetSaver<MODEL>::doFinalize(const State_ & xx) {
  Log::info() << "StateSetSaver::doFinalize (empty) done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops
