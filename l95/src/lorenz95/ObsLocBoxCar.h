/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSLOCBOXCAR_H_
#define LORENZ95_OBSLOCBOXCAR_H_

#include <ostream>

#include "oops/base/ObsLocalizationBase.h"

#include "lorenz95/L95Traits.h"
#include "lorenz95/ObsLocParameters.h"

namespace lorenz95 {
// Forward declarations
  class Iterator;
  class ObsTable;
  class ObsVec1D;

/// Observation space localization for BoxCar function
class ObsLocBoxCar: public oops::ObsLocalizationBase<L95Traits, L95ObsTraits> {
 public:
  typedef ObsLocParameters Parameters_;

  ObsLocBoxCar(const Parameters_ &, const ObsTable &);

  /// compute localization and update localization values in \p locfactor
  /// (missing value is for obs outside of localization)
  void computeLocalization(const Iterator &, ObsVec1D & locfactor) const override;

 private:
  void print(std::ostream &) const override;

  /// Gaspari-Cohn localization distance (localization goes to zero at rscale_)
  const double rscale_;

  /// ObsSpace associated with the observations
  const ObsTable & obsdb_;
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSLOCBOXCAR_H_
