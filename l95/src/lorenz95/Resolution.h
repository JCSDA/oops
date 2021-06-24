/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_RESOLUTION_H_
#define LORENZ95_RESOLUTION_H_

#include <iostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "lorenz95/Iterator.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace lorenz95 {

class Iterator;

// -----------------------------------------------------------------------------
/// \brief Parameters controlling a Lorenz95 model's resolution.
class ResolutionParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ResolutionParameters, Parameters)

 public:
  /// \brief Number of gridpoints.
  oops::RequiredParameter<int> resol{"resol", this};
};

// -----------------------------------------------------------------------------
/// Handles resolution.

class Resolution : public util::Printable {
 public:
  typedef ResolutionParameters Parameters_;

  Resolution(const ResolutionParameters & parameters, const eckit::mpi::Comm & comm)
              : resol_(parameters.resol), comm_(comm)
    {ASSERT(comm_.size() == 1);}
  explicit Resolution(const int resol): resol_(resol), comm_(oops::mpi::myself())
    {ASSERT(comm_.size() == 1);}

  int npoints() const {return resol_;}

  Iterator begin() const;
  Iterator end() const;
  std::vector<double> verticalCoord(std::string &) const;
  std::vector<size_t> variableSizes(const oops::Variables &) const;
  const eckit::mpi::Comm & getComm() const {return comm_;}

 private:
  void print(std::ostream & os) const {os << resol_;}
  const int resol_;
  const eckit::mpi::Comm & comm_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_RESOLUTION_H_
