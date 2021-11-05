/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/ObsErrorDiag.h"

namespace oops {

void ObsErrorDiagParameters::deserialize(util::CompositePath &path,
                                         const eckit::Configuration &config) {
  ObsErrorParametersBase::deserialize(path, config);

  if (zeroMeanPerturbations) {
    if (member.value() == boost::none)
      throw eckit::UserError(
          path.path() +
          ": The 'member' option must be set when 'zero-mean perturbations' is set to true",
          Here());
    if (numberOfMembers.value() == boost::none)
      throw eckit::UserError(
          path.path() +
          ": The 'number of members' option must be set when 'zero-mean perturbations' is set "
          "to true", Here());
    if (*member.value() > *numberOfMembers.value())
      throw eckit::UserError(
          path.path() +
          ": 'member' must not be greater than 'number of members'", Here());
  }
}

}  // namespace oops
