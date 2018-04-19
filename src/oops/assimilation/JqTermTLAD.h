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

 public:
  explicit JqTermTLAD(unsigned nsub) : mxi_(), nsubwin_(nsub) {}
  explicit JqTermTLAD(unsigned, const Increment4D_ &);
  ~JqTermTLAD() {}

  void clear() {xi_.clear();}
  void computeModelErrorTL(Increment4D_ &);

  GeneralizedDepartures * releaseOutputFromTL() override {return 0;}

 private:
  void doInitializeTraj(const State_ &, const util::DateTime &,
                      const util::Duration &) override {}
  void doProcessingTraj(const State_ &) override {}
  void doFinalizeTraj(const State_ &) override {}

  void doInitializeTL(const Increment_ &, const util::DateTime &,
                      const util::Duration &) override {}
  void doProcessingTL(const Increment_ &) override {}
  void doFinalizeTL(const Increment_ &) override;

  void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(Increment_ &) override {}
  void doLastAD(Increment_ &) override {}

  std::vector<Increment_> mxi_;
  std::vector<Increment_> xi_;
  const unsigned nsubwin_;
  unsigned current_;
};

// =============================================================================

template <typename MODEL>
JqTermTLAD<MODEL>::JqTermTLAD(unsigned nsub, const Increment4D_ & dx)
  : mxi_(), xi_(), nsubwin_(nsub), current_(0)
{
  for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
    xi_.push_back(dx[jsub]);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFinalizeTL(const Increment_ & dx) {
  if (mxi_.size() < nsubwin_ - 1) mxi_.push_back(dx);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::computeModelErrorTL(Increment4D_ & dx) {
// Compute x_i - M(x_{i-1})
  for (unsigned jsub = 1; jsub < nsubwin_; ++jsub) {
    dx[jsub] -= mxi_[jsub-1];
    Log::info() << "CostJbJq: x_" << jsub << " - M(x_" << jsub-1 << ")" << dx[jsub] << std::endl;
  }
  mxi_.clear();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTLAD<MODEL>::doFirstAD(Increment_ & dx, const util::DateTime &,
                                const util::Duration &) {
  if (current_ > 0) {
    dx -= xi_.back();
    xi_.pop_back();
  }
  current_ += 1;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_JQTERMTLAD_H_
