/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSLOCGC99_H_
#define LORENZ95_OBSLOCGC99_H_

#include <ostream>

#include "eckit/config/Configuration.h"
#include "oops/util/Printable.h"

#include "lorenz95/ObsData1D.h"

// Forward declarations
namespace lorenz95 {
  class Iterator;
  class ObsTable;
  class ObsVec1D;

/// Observation space localization for Lorenz 95 model (Gaspari-Cohn)
class ObsLocGC99: public util::Printable {
 public:
  ObsLocGC99(const eckit::Configuration &, const ObsTable &);

  /// compute localization and save localization values in \p obsvector and
  /// localization flags (1: outside of localization; 0: inside localization area)
  /// in \p outside
  void computeLocalization(const Iterator &, ObsData1D<int> & outside,
                           ObsVec1D & obsvector) const;

 private:
  void print(std::ostream &) const override;

  /// Gaspari-Cohn localization distance (localization goes to zero at rscale_)
  const double rscale_;

  /// ObsSpace associated with the observations
  const ObsTable & obsdb_;
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSLOCGC99_H_
