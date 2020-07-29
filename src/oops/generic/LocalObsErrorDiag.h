/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_LOCALOBSERRORDIAG_H_
#define OOPS_GENERIC_LOCALOBSERRORDIAG_H_

#include <memory>
#include <sstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/base/ObsErrorBase.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/generic/ObsErrorDiag.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Diagonal observation error covariance matrix with R-localization
//  localization (inflating observation error variances) is done in the constructor
//  the rest of the methods are not overriden; ObsErrorDiag methods would be used
//  instead

template<typename OBS>
class LocalObsErrorDiag : public ObsErrorDiag<OBS> {
  typedef ObsLocalizationBase<OBS>   ObsLocalization_;
  typedef ObsSpace<OBS>              ObsSpace_;

 public:
/// Initialize and inflate local R for obs. localization
  LocalObsErrorDiag(const eckit::Configuration &, const ObsSpace_ &);

 private:
  void print(std::ostream &) const override;
};

// =============================================================================

template<typename OBS>
LocalObsErrorDiag<OBS>::LocalObsErrorDiag
    (const eckit::Configuration & conf, const ObsSpace_ & obsdb)
  : ObsErrorDiag<OBS>(conf, obsdb)
{
// if Localization section is available; localize covariances
  if (conf.has("localization")) {
    eckit::LocalConfiguration locconf(conf, "localization");
    std::unique_ptr<ObsLocalization_> loc(ObsLocalizationFactory<OBS>::create(locconf, obsdb));
    loc->multiply(this->inverseVariance_);
  }
}

// -----------------------------------------------------------------------------

template<typename OBS>
void LocalObsErrorDiag<OBS>::print(std::ostream & os) const {
  os << "Localized diagonal observation error covariance, inverse variances: "
     << this->inverseVariance_ << std::endl;
}


// -----------------------------------------------------------------------------


}  // namespace oops

#endif  // OOPS_GENERIC_LOCALOBSERRORDIAG_H_
