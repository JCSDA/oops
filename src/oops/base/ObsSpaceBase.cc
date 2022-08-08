/*
 * (C) Copyright 2018  UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/ObsSpaceBase.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
int ObsSpaceBase::instances_ = 0;
// -----------------------------------------------------------------------------

ObsSpaceBase::ObsSpaceBase(const ObsSpaceParametersBase & params, const eckit::mpi::Comm & comm,
                           const util::DateTime & bgn, const util::DateTime & end)
  : winbgn_(bgn), winend_(end), instance_(++instances_) {
//
// Determine seed for random number generator that is reproducible when re-running
// but does not repeat itself over analysis cycles, ensemble members or obs type
//
  util::DateTime ref(1623, 6, 19, 0, 0, 0);
  ASSERT(winbgn_ > ref);
  util::Duration dt(winbgn_ - ref);
  seed_ = dt.toSeconds();

  // Won't repeat if more seconds between analysis cycles than members in EDA
  seed_ += params.obsPerturbationsSeed;

//             31622400 seconds max in 1 year
//          12197962800 seed at this step for 2010-01-01T03:00:00Z
  seed_ += 100000000000 * instance_;
}

// -----------------------------------------------------------------------------

}  // namespace oops
