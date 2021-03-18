/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_MINIMIZERUTILS_H_
#define OOPS_ASSIMILATION_MINIMIZERUTILS_H_

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/assimilation/ControlIncrement.h"

namespace oops {

/// Prints to Log::info gradient reduction \p grad and normalized gradient reduction \p norm
/// for iteration \p iteration
void printNormReduction(int iteration, const double & grad, const double & norm);
/// Prints to Log::info cost function values for \p costJ, \p costJb, \p costJoJc for
/// iteration \p iteration
void printQuadraticCostFunction(int iteration, const double & costJ, const double & costJb,
                                const double & costJoJc);

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void writeIncrement(const eckit::Configuration & config,
                    const ControlIncrement<MODEL, OBS> & dx, const int & loop) {
// Write out the increment

  if (config.has("online diagnostics")) {
    const eckit::LocalConfiguration onlineDiag(config, "online diagnostics");
    bool writeinc = onlineDiag.getBool("write increment", false);

    if (writeinc) {
      // print log
      Log::info() << "Write Increment - starting: " << loop << std::endl << std::endl;

      const eckit::LocalConfiguration incConf(config, "increment");

      // write increment
      dx.write(incConf);

      // print log
      Log::info() << std::endl << "Write Increment: done." << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void writeKrylovBasis(const eckit::Configuration & config,
                        const ControlIncrement<MODEL, OBS> & dx,
                        const int & loop) {
// Write out the increment
  if (config.has("online diagnostics")) {
    eckit::LocalConfiguration diagConf(config, "online diagnostics");
    bool writeinc = diagConf.getBool("write basis", false);

    if (writeinc) {
      // print log
      Log::info() << "Write Krylov Basis: starting: " << loop << std::endl;

      eckit::LocalConfiguration basisConf(config, "krylov basis");
      basisConf.set("iteration", loop);

      // write increment
      dx.write(basisConf);

      // print log
      Log::info() << "Write Krylov Basis: done." << std::endl;
    }
  }
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_MINIMIZERUTILS_H_
