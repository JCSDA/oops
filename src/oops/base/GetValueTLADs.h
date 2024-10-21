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

#include "oops/base/GetValuePosts.h"
#include "oops/base/GetValues.h"
#include "oops/base/Increment.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/base/State.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/VariableChange.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Computes observation equivalent TL and AD to/from increments.

template <typename MODEL, typename OBS>
class GetValueTLADs : public PostBaseTLAD<MODEL> {
  typedef Increment<MODEL>          Increment_;
  typedef State<MODEL>              State_;
  typedef VariableChange<MODEL>     VariableChange_;
  typedef typename VariableChange_::Parameters_ VariableChangeParameters_;
  typedef std::shared_ptr<GetValues<MODEL, OBS>> GetValPtr_;
  typedef std::unique_ptr<LinearVariableChange<MODEL>> CVarPtr_;

 public:
  GetValueTLADs(const GetValuesParameters<MODEL> &,
                const util::DateTime &, const util::DateTime &);

  void append(GetValPtr_);

  GetValPtr_ & operator[](const std::size_t ii) {return getvals_[ii];}

 private:
// Methods
  void doInitializeTraj(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingTraj(const State_ &) override;
  void doFinalizeTraj(const State_ &) override;

  void doInitializeTL(const Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingTL(const Increment_ &) override;
  void doFinalizeTL(const Increment_ &) override;

  void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(Increment_ &) override;
  void doLastAD(Increment_ &) override;

// Data
  Variables geovars_;
  Variables linvars_;
  std::vector<GetValPtr_> getvals_;
  std::map<util::DateTime, CVarPtr_> chvartlad_;
  GetValuesParameters<MODEL> params_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
GetValueTLADs<MODEL, OBS>::GetValueTLADs(const GetValuesParameters<MODEL> & params,
              const util::DateTime & bgn, const util::DateTime & end)
  : PostBaseTLAD<MODEL>(bgn, end), geovars_(), linvars_(), getvals_(), chvartlad_(),
    params_(params)
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
  for (GetValPtr_ getval : getvals_) getval->initialize(tstep);
  Log::trace() << "GetValueTLADs::doInitializeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doProcessingTraj(const State_ & xx) {
  Log::trace() << "GetValueTLADs::doProcessingTraj start" << std::endl;

  VariableChange_ chvar(params_.variableChange.value(), xx.geometry());

  State_ zz(xx);
  chvar.changeVar(zz, geovars_);

  // Optimization: if underlying interpolations are done using the atlas representation of the
  // model data, then call FieldSet::haloExchange before entering the loop over obs types:
  PreProcessHelper<MODEL>::preProcessModelData(zz);

  for (GetValPtr_ getval : getvals_) getval->process(zz);

  CVarPtr_ cvtlad(new LinearVariableChange<MODEL>(xx.geometry(),
                  params_.variableChange.value().toConfiguration()));
  cvtlad->changeVarTraj(xx, linvars_);
  chvartlad_[xx.validTime()] = std::move(cvtlad);

  Log::trace() << "GetValueTLADs::doProcessingTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doFinalizeTraj(const State_ &) {
  Log::trace() << "GetValueTLADs::doFinalizeTraj start" << std::endl;
  for (GetValPtr_ getval : getvals_) getval->finalize();
  Log::trace() << "GetValueTLADs::doFinalizeTraj done" << std::endl;
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

  Increment_ dz(dx);
  chvartlad_[now]->changeVarTL(dz, linvars_);

  // Optimization: if underlying interpolations are done using the atlas representation of the
  // model data, then call FieldSet::haloExchange before entering the loop over obs types:
  PreProcessHelper<MODEL>::preProcessModelData(dz);

  for (GetValPtr_ getval : getvals_) getval->processTL(dz);

  Log::trace() << "GetValueTLADs::doProcessingTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doFinalizeTL(const Increment_ &) {
  Log::trace() << "GetValueTLADs::doFinalizeTL start" << std::endl;
  for (GetValPtr_ getval : getvals_) getval->finalizeTL();
  Log::trace() << "GetValueTLADs::doFinalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doFirstAD(Increment_ &, const util::DateTime &,
                                          const util::Duration & tstep) {
  Log::trace() << "GetValueTLADs::doFirstAD start" << std::endl;
  for (GetValPtr_ getval : getvals_) getval->finalizeAD(tstep);
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

  for (GetValPtr_ getval : getvals_) {
    getval->processAD(dz);
  }

  // Optimization: if underlying interpolations are done using the atlas representation of the
  // model data, then call FieldSet::adjointHaloExchange after exiting the loop over obs types:
  PreProcessHelper<MODEL>::preProcessModelDataAD(dz);

  dz.synchronizeFields();

  chvartlad_[now]->changeVarAD(dz, dx.variables());
  dx += dz;

  Log::trace() << "GetValueTLADs::doProcessingAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValueTLADs<MODEL, OBS>::doLastAD(Increment_ &) {
  Log::trace() << "GetValueTLADs::doLastAD start" << std::endl;
  for (GetValPtr_ getval : getvals_) getval->initializeAD();
  Log::trace() << "GetValueTLADs::doLastAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GETVALUETLADS_H_
