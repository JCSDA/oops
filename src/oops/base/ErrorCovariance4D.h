/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_ERRORCOVARIANCE4D_H_
#define OOPS_BASE_ERRORCOVARIANCE4D_H_

#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/ErrorCovariance.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic model space error 4D covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class ErrorCovariance4D : public ModelSpaceCovariance4DBase<MODEL> {
  typedef Geometry<MODEL>        Geometry_;
  typedef ErrorCovariance<MODEL> ErrorCovariance_;
  typedef Increment<MODEL>       Increment_;
  typedef Increment4D<MODEL>     Increment4D_;
  typedef State4D<MODEL>         State4D_;

 public:
  ErrorCovariance4D(const Geometry_ &, const Variables &,
                    const eckit::Configuration &, const State4D_ &, const State4D_ &);
  ~ErrorCovariance4D();

  void randomize(Increment4D_ &) const;

 private:
  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  unsigned nsubwin_;
  boost::ptr_vector<ErrorCovariance_> cov_;
};

// =============================================================================

/// Constructor, destructor
// -----------------------------------------------------------------------------
template<typename MODEL>
ErrorCovariance4D<MODEL>::ErrorCovariance4D(const Geometry_ & resol, const Variables & vars,
                                            const eckit::Configuration & conf,
                                            const State4D_ & xb, const State4D_ & fg)
  : ModelSpaceCovariance4DBase<MODEL>(xb, fg, resol, conf),
    cov_()
{
  Log::trace() << "ErrorCovariance4D::ErrorCovariance4D start" << std::endl;
// No cross-timeslot covariances
  std::vector < eckit::LocalConfiguration > confTime;
  conf.get("covariance_time", confTime);
  nsubwin_ = confTime.size();
  ASSERT(xb.size() == nsubwin_);
  for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
    cov_.push_back(new ErrorCovariance_(resol, vars, confTime[jsub], xb[jsub], fg[jsub]));
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ErrorCovariance4D<MODEL>::~ErrorCovariance4D() {
  Log::trace() << "ErrorCovariance4D destructed." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ErrorCovariance4D<MODEL>::doRandomize(Increment4D_ & dx) const {
// No cross-timeslot covariances
  for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
    int isub = jsub+dx.first();
    cov_[jsub].randomize(dx[isub]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ErrorCovariance4D<MODEL>::doMultiply(const Increment4D_ & dxi,
                                          Increment4D_ & dxo) const {
// No cross-timeslot covariances
  for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
    int isub = jsub+dxi.first();
    cov_[jsub].multiply(dxi[isub], dxo[isub]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ErrorCovariance4D<MODEL>::doInverseMultiply(const Increment4D_ & dxi,
                                                 Increment4D_ & dxo) const {
  IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_ERRORCOVARIANCE4D_H_
