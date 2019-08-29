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

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

#include "lorenz95/ObsTableView.h"

namespace lorenz95 {
  class LocsL95;

// -----------------------------------------------------------------------------

class ObsDiags1D : public util::Printable {
 public:
  ObsDiags1D(const ObsTableView &, const LocsL95 &, const oops::Variables &) {}
  ~ObsDiags1D() {}

// I/O
  void save(const std::string &) const {}

 private:
  void print(std::ostream &) const {}
};
//// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSDIAGS1D_H_
