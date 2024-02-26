/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2023 UCAR
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

#include "oops/assimilation/JqTerm.h"
#include "oops/base/Increment.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/base/State.h"
#include "oops/mpi/mpi.h"

namespace util {
class DateTime;
class Duration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class JqTermTLAD : public PostBaseTLAD<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;
  typedef JqTerm<MODEL>              JqTerm_;

 public:
  explicit JqTermTLAD(const eckit::mpi::Comm &);
  ~JqTermTLAD() {}

  std::shared_ptr<JqTerm_> getJq() {return jq_;}
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
  std::shared_ptr<JqTerm_> jq_;
};

// =============================================================================

template <typename MODEL>
JqTermTLAD<MODEL>::JqTermTLAD(const eckit::mpi::Comm & comm)
  : commTime_(comm), jq_(new JqTerm_())
{
  Log::trace() << "JqTermTLAD::JqTermTLAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFinalizeTraj(const State_ & xx) {
  Log::trace() << "JqTermTLAD::doFinalizeTraj start" << std::endl;
  jq_->finalize(xx);
  Log::trace() << "JqTermTLAD::doFinalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFinalizeTL(const Increment_ & dx) {
  Log::trace() << "JqTermTLAD::doFinalizeTL start" << std::endl;
  size_t mytime = commTime_.rank();
  if (mytime + 1 < commTime_.size()) oops::mpi::send(commTime_, dx, mytime+1, 2468);
  Log::trace() << "JqTermTLAD::doFinalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::computeModelErrorTL(Increment_ & dx) {
  Log::trace() << "JqTermTLAD::computeModelErrorTL start" << std::endl;
// Compute x_i - M(x_{i-1})
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
