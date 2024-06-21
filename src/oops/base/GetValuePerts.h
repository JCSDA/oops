/*
 * (C) Copyright 2024 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_GETVALUEPERTS_H_
#define OOPS_BASE_GETVALUEPERTS_H_

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/GetValuePosts.h"
#include "oops/base/GetValueTLADs.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {


/// \brief Fills GeoVaLs with requested variables at requested locations during model run
template <typename MODEL, typename OBS>
class GetValuePerts : public GetValuePosts<MODEL, OBS> {
  typedef State<MODEL>                           State_;
  typedef GetValuePosts<MODEL, OBS>              GetValuePosts_;
  typedef GetValueTLADs<MODEL, OBS>              GetValueTLADs_;
  typedef Increment<MODEL>                       Increment_;
  typedef LinearVariableChange<MODEL>            LinearVariableChange_;
  typedef VariableChange<MODEL>     VariableChange_;
  typedef typename VariableChange_::Parameters_ VariableChangeParameters_;
  typedef std::shared_ptr<GetValues<MODEL, OBS>> GetValuePtr_;

 public:
/// \brief Saves Locations and Variables to be processed
  GetValuePerts(const GetValuesParameters<MODEL> &, const std::shared_ptr<GetValueTLADs_> &,
                const Variables &);

  void append(GetValuePtr_) override {};
  void clear() override {};

 private:
/// \brief called at each model step: fill in GeoVaLs for the current time slot
  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;

  std::shared_ptr<GetValueTLADs_> getValTL_;
  Variables incVars_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValuePerts<MODEL, OBS>::GetValuePerts(const GetValuesParameters<MODEL>& params,
                                         const std::shared_ptr<GetValueTLADs_> & getValueTLADs,
                                         const Variables & incVars)
  : GetValuePosts_(params), getValTL_(), incVars_(incVars)
{
  getValTL_ = getValueTLADs;
  Log::trace() << "GetValuePerts::GetValuePerts" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePerts<MODEL, OBS>::doInitialize(const State_ & xx, const util::DateTime & dt,
                                             const util::Duration & tstep) {
  Log::trace() << "GetValuePerts::doInitialize start" << std::endl;
  Increment_ zz_inc(xx.geometry(), incVars_, xx.validTime());
  getValTL_->initializeTL(zz_inc, dt, tstep);
  Log::trace() << "GetValuePerts::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePerts<MODEL, OBS>::doProcessing(const State_ & xx) {
  Log::trace() << "GetValuePerts::doProcessing start" << std::endl;

  State_ zz_zero(xx);
  zz_zero.zero();
  Increment_ zz_inc(xx.geometry(), incVars_, xx.validTime());
  zz_inc.diff(xx, zz_zero);

  getValTL_->processTL(zz_inc);

  Log::trace() << "GetValuePerts::doProcessing done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePerts<MODEL, OBS>::doFinalize(const State_ & xx) {
  Log::trace() << "GetValuePosts::doFinalize start" << std::endl;
  Increment_ zz_inc(xx.geometry(), incVars_, xx.validTime());
  getValTL_->finalizeTL(zz_inc);
  Log::trace() << "GetValuePosts::doFinalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GETVALUEPERTS_H_
