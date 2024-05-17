/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_OBSDIAGS1D_H_
#define LORENZ95_OBSDIAGS1D_H_

#include <ostream>
#include <string>

#include "oops/base/ObsVariables.h"
#include "oops/util/Printable.h"

#include "lorenz95/ObsTable.h"

namespace oops {
  template <typename OBS> class Locations;
}

namespace lorenz95 {
  struct L95ObsTraits;

// -----------------------------------------------------------------------------

class ObsDiags1D : public util::Printable {
 public:
  ObsDiags1D(const ObsTable &, const oops::Locations<L95ObsTraits> &,
             const oops::ObsVariables &) {}
  ~ObsDiags1D() {}

// I/O
  void save(const std::string &) const {}

 private:
  void print(std::ostream &) const {}
};
//// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSDIAGS1D_H_
