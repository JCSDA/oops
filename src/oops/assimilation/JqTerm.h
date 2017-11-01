/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_JQTERM_H_
#define OOPS_ASSIMILATION_JQTERM_H_

#include <vector>

#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/PostBase.h"
#include "oops/interface/State.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class JqTerm : public PostBase< State<MODEL> > {
  typedef State<MODEL>          State_;
  typedef Increment4D<MODEL>    Increment4D_;
  typedef State4D<MODEL>        State4D_;

 public:
  explicit JqTerm(unsigned nsub) : mxi_(), nsubwin_(nsub) {}
  ~JqTerm() {}

  void linearize() {}
  void computeModelError(const State4D_ &, Increment4D_ &);

 private:
  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override {}
  void doProcessing(const State_ &) override {}
  void doFinalize(const State_ &) override;

  std::vector<State_> mxi_;
  const unsigned nsubwin_;
};

// =============================================================================

template <typename MODEL>
void JqTerm<MODEL>::doFinalize(const State_ & xx) {
  if (mxi_.size() < nsubwin_ - 1) mxi_.push_back(xx);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTerm<MODEL>::computeModelError(const State4D_ & fg, Increment4D_ & dx) {
// Compute x_i - M(x_{i-1})
  for (unsigned jsub = 1; jsub < nsubwin_; ++jsub) {
    dx[jsub].diff(fg[jsub], mxi_[jsub-1]);
    Log::info() << "CostJbJq: x_" << jsub << " - M(x_" << jsub-1 << ")" << dx[jsub] << std::endl;
  }
  mxi_.clear();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_JQTERM_H_
