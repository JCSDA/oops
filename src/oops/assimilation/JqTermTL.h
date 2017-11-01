/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_JQTERMTL_H_
#define OOPS_ASSIMILATION_JQTERMTL_H_

#include <vector>

#include "oops/base/PostBaseTL.h"
#include "oops/interface/Increment.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {
  template<typename MODEL> class Increment4D;

// -----------------------------------------------------------------------------

template <typename MODEL> class JqTermTL : public PostBaseTL<Increment<MODEL> > {
  typedef Increment<MODEL>           Increment_;
  typedef Increment4D<MODEL>         Increment4D_;

 public:
  explicit JqTermTL(unsigned nsub) : mxi_(), nsubwin_(nsub) {}
  ~JqTermTL() {}

  void computeModelErrorTL(Increment4D_ &);

  GeneralizedDepartures * releaseOutputFromTL() override {return 0;}

 private:
  void doInitializeTL(const Increment_ &, const util::DateTime &,
                      const util::Duration &) override {}
  void doProcessingTL(const Increment_ &) override {}
  void doFinalizeTL(const Increment_ &) override;

  std::vector<Increment_> mxi_;
  const unsigned nsubwin_;
};

// =============================================================================

template <typename MODEL>
void JqTermTL<MODEL>::doFinalizeTL(const Increment_ & dx) {
  if (mxi_.size() < nsubwin_ - 1) mxi_.push_back(dx);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermTL<MODEL>::computeModelErrorTL(Increment4D_ & dx) {
// Compute x_i - M(x_{i-1})
  for (unsigned jsub = 1; jsub < nsubwin_; ++jsub) {
    dx[jsub] -= mxi_[jsub-1];
    Log::info() << "CostJbJq: x_" << jsub << " - M(x_" << jsub-1 << ")" << dx[jsub] << std::endl;
  }
  mxi_.clear();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_JQTERMTL_H_
