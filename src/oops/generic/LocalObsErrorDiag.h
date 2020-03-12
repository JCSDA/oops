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

template<typename MODEL>
class LocalObsErrorDiag : public ObsErrorDiag<MODEL> {
  typedef ObsLocalizationBase<MODEL>   ObsLocalization_;
  typedef ObsSpace<MODEL>              ObsSpace_;

 public:
/// Initialize and inflate local R for obs. localization
  LocalObsErrorDiag(const eckit::Configuration &, const ObsSpace_ &);

 private:
  void print(std::ostream &) const override;
};

// =============================================================================

template<typename MODEL>
LocalObsErrorDiag<MODEL>::LocalObsErrorDiag
    (const eckit::Configuration & conf, const ObsSpace_ & obsdb)
  : ObsErrorDiag<MODEL>(conf, obsdb)
{
// if Localization section is available; localize covariances
  if (conf.has("Localization")) {
    eckit::LocalConfiguration locconf(conf, "Localization");
    std::unique_ptr<ObsLocalization_> loc(ObsLocalizationFactory<MODEL>::create(locconf, obsdb));
    loc->multiply(this->inverseVariance_);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalObsErrorDiag<MODEL>::print(std::ostream & os) const {
  os << "LocalObsErrorDiag<MODEL>::print not implemeted yet";
}


// -----------------------------------------------------------------------------


}  // namespace oops

#endif  // OOPS_GENERIC_LOCALOBSERRORDIAG_H_
