/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_JQTERMTLAD_H_
#define OOPS_ASSIMILATION_JQTERMTLAD_H_

#include <memory>
#include <vector>

#include "oops/base/PostBaseTLAD.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class JqTermTLAD : public PostBaseTLAD<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  explicit JqTermTLAD(const eckit::mpi::Comm &);
  ~JqTermTLAD() {}

  void clear() {xi_.reset();}
// void computeModelErrorTraj(const State_ &, Increment_ &);  // not used
  State_ & getMxi() const;
  void computeModelErrorTL(Increment_ &);

  void setupAD(const Increment_ & dx);

 private:
  void doInitializeTraj(const State_ &, const util::DateTime &,
                        const util::Duration &) override {}
  void doProcessingTraj(const State_ &) override {}
  void doFinalizeTraj(const State_ &) override;

  void doInitializeTL(const Increment_ &, const util::DateTime &,
                      const util::Duration &) override {}
  void doProcessingTL(const Increment_ &) override {}
  void doFinalizeTL(const Increment_ &) override;

  void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(Increment_ &) override {}
  void doLastAD(Increment_ &) override {}

  const eckit::mpi::Comm & commTime_;
  std::unique_ptr<State_> xtraj_;
  std::unique_ptr<Increment_> mxi_;
  std::unique_ptr<Increment_> xi_;
};

// =============================================================================

template <typename MODEL>
JqTermTLAD<MODEL>::JqTermTLAD(const eckit::mpi::Comm & comm)
  : commTime_(comm), xtraj_(), mxi_(), xi_()
{
  Log::trace() << "JqTermTLAD::JqTermTLAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFinalizeTraj(const State_ & xx) {
  Log::trace() << "JqTermTLAD::doFinalizeTraj start" << std::endl;
  xtraj_.reset(new State_(xx));
  Log::trace() << "JqTermTLAD::doFinalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

/*
template <typename MODEL>
void JqTermTLAD<MODEL>::computeModelErrorTraj(const State_ & fg, Increment_ & dx) {
  Log::trace() << "JqTermTLAD::computeModelErrorTraj start" << std::endl;

  static int tag = 83655;
  size_t mytime = commTime_.rank();
// Send values of M(x_i) at end of my subwindow to next subwindow
  if (mytime + 1 < commTime_.size()) {
    Log::debug() << "JqTermTLAD::computeModelErrorTraj: sending to " << mytime+1
                 << " " << tag << std::endl;
    oops::mpi::send(commTime_, fg, mytime+1, tag);
    Log::debug() << "JqTermTLAD::computeModelErrorTraj: sent to " << mytime+1
                 << " " << tag << std::endl;
  }

// Receive values at beginning of my subwindow from previous subwindow
  if (mytime > 0) {
    State_ mxi(fg);
    Log::debug() << "JqTermTLAD::computeModelErrorTraj: receiving from " << mytime-1
                 << " " << tag << std::endl;
    oops::mpi::receive(commTime_, mxi, mytime-1, tag);
    Log::debug() << "JqTermTLAD::computeModelErrorTraj: received from " << mytime-1
                 << " " << tag << std::endl;

//  Compute x_i - M(x_{i-1})
    dx.diff(fg, mxi);
  }
  ++tag;
  Log::trace() << "JqTermTLAD::computeModelErrorTraj done" << std::endl;
}
*/

// -----------------------------------------------------------------------------

template <typename MODEL>
State<MODEL> & JqTermTLAD<MODEL>::getMxi() const {
  Log::trace() << "JqTermTLAD::getMxi" << std::endl;
// Retrieve M(x-i)
  return *xtraj_;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFinalizeTL(const Increment_ & dx) {
  Log::trace() << "JqTermTLAD::doFinalizeTL start" << std::endl;
  Log::debug() << "JqTermTLAD::doFinalizeTL MPI size " << commTime_.size() << std::endl;
  Log::debug() << "JqTermTLAD::doFinalizeTL MPI rank " << commTime_.rank() << std::endl;
  size_t mytime = commTime_.rank();
  if (mytime + 1 < commTime_.size()) oops::mpi::send(commTime_, dx, mytime+1, 2468);
  Log::trace() << "JqTermTLAD::doFinalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::computeModelErrorTL(Increment_ & dx) {
  Log::trace() << "JqTermTLAD::computeModelErrorTL start" << std::endl;
// Compute x_i - M(x_{i-1})
  Log::debug() << "JqTermTLAD::computeModelErrorTL MPI size " << commTime_.size() << std::endl;
  Log::debug() << "JqTermTLAD::computeModelErrorTL MPI rank " << commTime_.rank() << std::endl;
  size_t mytime = commTime_.rank();
  if (mytime > 0) {
    Increment_ mxim1(dx, false);
    oops::mpi::receive(commTime_, mxim1, mytime-1, 2468);
    dx -= mxim1;
  }
  Log::info() << "JqTermTLAD: x_i - M(x_i)" << dx << std::endl;
  Log::trace() << "JqTermTLAD::computeModelErrorTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::setupAD(const Increment_ & dx) {
  Log::trace() << "JqTermTLAD::setupAD start" << std::endl;
  size_t mytime = commTime_.rank();
  if (mytime > 0) oops::mpi::send(commTime_, dx, mytime-1, 8642);
  Log::trace() << "JqTermTLAD::setupAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFirstAD(Increment_ & dx, const util::DateTime &,
                                  const util::Duration &) {
  Log::trace() << "JqTermTLAD::doFirstAD start" << std::endl;
  size_t mytime = commTime_.rank();
  if (mytime + 1 < commTime_.size()) {
    Increment_ xip1(dx, false);
    oops::mpi::receive(commTime_, xip1, mytime+1, 8642);
    dx -= xip1;
  }
  Log::trace() << "JqTermTLAD::doFirstAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_JQTERMTLAD_H_
