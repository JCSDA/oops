/*
 * (C) Copyright 2021-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_GETVALUEPOSTS_H_
#define OOPS_BASE_GETVALUEPOSTS_H_

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/GetValues.h"
#include "oops/base/PostBase.h"
#include "oops/base/State.h"
#include "oops/interface/VariableChange.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// \brief Fills GeoVaLs with requested variables at requested locations during model run
template <typename MODEL, typename OBS>
class GetValuePosts : public PostBase<State<MODEL>> {
  typedef VariableChange<MODEL>     VariableChange_;
  typedef typename VariableChange_::Parameters_ VariableChangeParameters_;
  typedef State<MODEL>              State_;
  typedef std::shared_ptr<GetValues<MODEL, OBS>> GetValuePtr_;

 public:
/// \brief Saves Locations and Variables to be processed
  GetValuePosts();

  void append(GetValuePtr_);

 private:
/// \brief initialization before model run: sets up GetValues and allocate GeoVaLs
  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
/// \brief called at each model step: fill in GeoVaLs for the current time slot
  void doProcessing(const State_ &) override;

// Data
  std::vector<GetValuePtr_> getvals_;
  Variables geovars_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValuePosts<MODEL, OBS>::GetValuePosts() : PostBase<State_>(), getvals_(), geovars_() {
  Log::trace() << "GetValuePosts::GetValuePosts" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePosts<MODEL, OBS>::append(GetValuePtr_ getval) {
  Log::trace() << "GetValuePosts::append start" << std::endl;
  getvals_.push_back(getval);
  geovars_ += getval->requiredVariables();
  Log::trace() << "GetValuePosts::append done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePosts<MODEL, OBS>::doInitialize(const State_ &, const util::DateTime &,
                                             const util::Duration & tstep) {
  Log::trace() << "GetValuePosts::doInitialize start" << std::endl;
  for (GetValuePtr_ getval : getvals_) getval->initialize(tstep);
  Log::trace() << "GetValuePosts::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePosts<MODEL, OBS>::doProcessing(const State_ & xx) {
  Log::trace() << "GetValuePosts::doProcessing start" << std::endl;

  eckit::LocalConfiguration chvarconf;
  VariableChangeParameters_ params;
  params.validateAndDeserialize(chvarconf);
  VariableChange_ chvar(params, xx.geometry());

  State_ zz(xx);
  chvar.changeVar(zz, geovars_);

  for (GetValuePtr_ getval : getvals_) getval->process(zz);

  Log::trace() << "GetValuePosts::doProcessing done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GETVALUEPOSTS_H_
