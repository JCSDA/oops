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

#include <vector>

#include "oops/assimilation/JqTerm.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {
  template<typename MODEL> class Increment4D;

// -----------------------------------------------------------------------------

template <typename MODEL>
class JqTermTLAD : public PostBaseTLAD<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef State<MODEL>               State_;
  typedef State4D<MODEL>        State4D_;

 public:
  explicit JqTermTLAD(unsigned nsub);
  ~JqTermTLAD() {}

  void clear() {xi_.clear();}
  void computeModelErrorTraj(const State4D_ &, Increment4D_ &);
  void computeModelErrorTL(Increment4D_ &);

  GeneralizedDepartures * releaseOutputFromTL() override {return 0;}
  void setupAD(const Increment4D_ & dx);

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

  JqTerm<MODEL> jq_;
  std::vector<Increment_> mxi_;
  std::vector<Increment_> xi_;
  const unsigned nsubwin_;
  unsigned current_;
};

// =============================================================================

template <typename MODEL>
JqTermTLAD<MODEL>::JqTermTLAD(unsigned nsub)
  : jq_(nsub), mxi_(), xi_(), nsubwin_(nsub), current_(0)
{
  Log::trace() << "JqTermTLAD::JqTermTLAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFinalizeTraj(const State_ & xx) {
  Log::trace() << "JqTermTLAD::doFinalizeTraj start" << std::endl;
  jq_.finalize(xx);
  Log::trace() << "JqTermTLAD::doFinalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::computeModelErrorTraj(const State4D_ & fg, Increment4D_ & dx) {
  Log::trace() << "JqTermTLAD::computeModelErrorTraj start" << std::endl;
  jq_.computeModelError(fg, dx);
  Log::trace() << "JqTermTLAD::computeModelErrorTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFinalizeTL(const Increment_ & dx) {
  Log::trace() << "JqTermTLAD::doFinalizeTL start" << std::endl;
  if (mxi_.size() < nsubwin_ - 1) mxi_.push_back(dx);
  Log::trace() << "JqTermTLAD::doFinalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::computeModelErrorTL(Increment4D_ & dx) {
  Log::trace() << "JqTermTLAD::computeModelErrorTL start" << std::endl;
// Compute x_i - M(x_{i-1})
  for (unsigned jsub = 1; jsub < nsubwin_; ++jsub) {
    dx[jsub] -= mxi_[jsub-1];
    Log::info() << "JqTermTLAD: x_" << jsub << " - M(x_" << jsub-1 << ")" << dx[jsub] << std::endl;
  }
  mxi_.clear();
  Log::trace() << "JqTermTLAD::computeModelErrorTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::setupAD(const Increment4D_ & dx) {
  Log::trace() << "JqTermTLAD::setupAD start" << std::endl;
  xi_.clear();
  for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
    xi_.push_back(dx[jsub]);
  }
  current_ = nsubwin_;
  Log::trace() << "JqTermTLAD::setupAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFirstAD(Increment_ & dx, const util::DateTime &,
                                  const util::Duration &) {
  Log::trace() << "JqTermTLAD::doFirstAD start" << std::endl;
  ASSERT(current_ >= 0);
  ASSERT(current_ <= nsubwin_);
  if (current_ < xi_.size()) {
    dx -= xi_[current_];
  }
  current_ -= 1;
  Log::trace() << "JqTermTLAD::doFirstAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_JQTERMTLAD_H_
