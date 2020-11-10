/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSLOCGC99_H_
#define LORENZ95_OBSLOCGC99_H_

#include <ostream>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

#include "lorenz95/L95Traits.h"

// Forward declarations
namespace lorenz95 {
  class ObsVec1D;

/// ObsLocalization matrix for Lorenz 95 model.

// -----------------------------------------------------------------------------
class ObsLocGC99: public util::Printable {
 public:
  static const std::string classname() {return "lorenz95::ObsLocGC99";}

  ObsLocGC99(const eckit::Configuration &, const ObsTableView &);
  void multiply(ObsVec1D &) const;

 private:
  void print(std::ostream &) const override;
  const ObsTableView & obsdb_;
  const double rscale_;
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSLOCGC99_H_
