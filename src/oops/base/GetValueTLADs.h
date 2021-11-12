/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_GETVALUETLADS_H_
#define OOPS_BASE_GETVALUETLADS_H_

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "oops/base/GetValueTLAD.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/base/State.h"
#include "oops/interface/ChangeVariables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.

template <typename MODEL, typename OBS>
class GetValueTLADs : public PostBaseTLAD<MODEL> {
  typedef ChangeVariables<MODEL>           ChangeVar_;
  typedef Increment<MODEL>                 Increment_;
  typedef State<MODEL>                     State_;
  typedef std::shared_ptr<GetValueTLAD<MODEL, OBS>> GetValPtr_;
  typedef std::unique_ptr<LinearVariableChangeBase<MODEL>> CVarPtr_;

 public:
  GetValueTLADs(const util::DateTime &, const util::DateTime &);

  void append(GetValPtr_);

 private:
// Methods
  void doInitializeTraj(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingTraj(const State_ &) override;
  void doFinalizeTraj(const State_ &) override {}

  void doInitializeTL(const Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingTL(const Increment_ &) override;
  void doFinalizeTL(const Increment_ &) override {}

  void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(Increment_ &) override;
  void doLastAD(Increment_ &) override {}

// Data
  Variables geovars_;
  Variables linvars_;
  std::vector<GetValPtr_> getvals_;
  std::map<util::DateTime, CVarPtr_> chvartlad_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
GetValueTLADs<MODEL, OBS>::GetValueTLADs(const util::DateTime & bgn, const util::DateTime & end)
  : PostBaseTLAD<MODEL>(bgn, end), geovars_(), linvars_(), getvals_(), chvartlad_()
{
  Log::trace() << "GetValueTLADs::GetValueTLADs" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::append(GetValPtr_ getval) {
  Log::trace() << "GetValuePosts::append start" << std::endl;
  getvals_.push_back(getval);
  geovars_ += getval->requiredVariables();
  linvars_ += getval->linearVariables();
  Log::trace() << "GetValuePosts::append done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doInitializeTraj(const State_ &, const util::DateTime &,
                                                 const util::Duration & tstep) {
  Log::trace() << "GetValueTLADs::doInitializeTraj start" << std::endl;
  for (GetValPtr_ getval : getvals_) getval->initializeTraj(tstep);
  Log::trace() << "GetValueTLADs::doInitializeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doProcessingTraj(const State_ & xx) {
  Log::trace() << "GetValueTLADs::doProcessingTraj start" << std::endl;

  eckit::LocalConfiguration chvarconf;
  chvarconf.set("variable change", "default");

  ChangeVar_ chvar(chvarconf, xx.geometry(), xx.variables(), geovars_);
  State_ zz(xx.geometry(), geovars_, xx.validTime());
  chvar.changeVar(xx, zz);

  for (GetValPtr_ getval : getvals_) getval->processTraj(zz);

  CVarPtr_ cvtlad(LinearVariableChangeFactory<MODEL>::create(xx, xx, xx.geometry(), chvarconf));
  chvartlad_[xx.validTime()] = std::move(cvtlad);

  Log::trace() << "GetValueTLADs::doProcessingTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doInitializeTL(const Increment_ &, const util::DateTime &,
                                               const util::Duration & tstep) {
  Log::trace() << "GetValueTLADs::doInitializeTL start" << std::endl;
  for (GetValPtr_ getval : getvals_) getval->initializeTL(tstep);
  Log::trace() << "GetValueTLADs::doInitializeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doProcessingTL(const Increment_ & dx) {
  Log::trace() << "GetValueTLADs::doProcessingTL start" << std::endl;
  const util::DateTime now = dx.validTime();
  ASSERT(chvartlad_.find(now) != chvartlad_.end());

  Increment_ dz(dx.geometry(), linvars_, dx.validTime());
  chvartlad_[now]->multiply(dx, dz);

  for (GetValPtr_ getval : getvals_) getval->processTL(dz);
  Log::trace() << "GetValueTLADs::doProcessingTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doFirstAD(Increment_ &, const util::DateTime &,
                                          const util::Duration & tstep) {
  Log::trace() << "GetValueTLADs::doFirstAD start" << std::endl;
  for (GetValPtr_ getval : getvals_) getval->initializeAD(tstep);
  Log::trace() << "GetValueTLADs::doFirstAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doProcessingAD(Increment_ & dx) {
  Log::trace() << "GetValueTLADs::doProcessingAD start" << std::endl;
  const util::DateTime now = dx.validTime();
  ASSERT(chvartlad_.find(now) != chvartlad_.end());

  Increment_ dz(dx.geometry(), linvars_, dx.validTime());
  dz.zero();

  Increment_ tmpz(dz);
  for (GetValPtr_ getval : getvals_) {
    tmpz.zero();
    getval->processAD(tmpz);
    dz += tmpz;
  }

  Increment_ tmpx(dx);
  tmpx.zero();
  chvartlad_[now]->multiplyAD(dz, tmpx);
  dx += tmpx;

  Log::trace() << "GetValueTLADs::doProcessingAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GETVALUETLADS_H_
