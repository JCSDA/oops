/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_JQTERMAD_H_
#define OOPS_ASSIMILATION_JQTERMAD_H_

#include <vector>

#include "oops/base/PostBaseAD.h"
#include "oops/interface/Increment.h"

namespace oops {
  class DateTime;
  class Duration;
  template<typename MODEL> class Increment4D;

// -----------------------------------------------------------------------------

template <typename MODEL> class JqTermAD : public PostBaseAD<Increment<MODEL> > {
  typedef Increment<MODEL>           Increment_;
  typedef Increment4D<MODEL>         Increment4D_;

 public:
  explicit JqTermAD(unsigned, const Increment4D_ &);
  ~JqTermAD() {}

  void clear() {xi_.clear();}

 private:
  void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(Increment_ &) override {}
  void doLastAD(Increment_ &) override {}

  const unsigned nsubwin_;
  unsigned current_;
  std::vector<Increment_> xi_;
};

// =============================================================================

template <typename MODEL>
JqTermAD<MODEL>::JqTermAD(unsigned nsub, const Increment4D_ & dx)
  : nsubwin_(nsub), current_(0), xi_()
{
  for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
    xi_.push_back(dx[jsub]);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTermAD<MODEL>::doFirstAD(Increment_ & dx, const util::DateTime &,
                                const util::Duration &) {
  if (current_ > 0) {
    dx -= xi_.back();
    xi_.pop_back();
  }
  current_ += 1;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_JQTERMAD_H_
