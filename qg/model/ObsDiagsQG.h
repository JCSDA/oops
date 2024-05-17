/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_OBSDIAGSQG_H_
#define QG_MODEL_OBSDIAGSQG_H_

#include <ostream>
#include <string>

#include "oops/base/ObsVariables.h"
#include "oops/util/Printable.h"

#include "oops/qg/ObsSpaceQG.h"

namespace oops {
  template <typename OBS> class Locations;
}

namespace qg {
  struct QgObsTraits;

// -----------------------------------------------------------------------------

class ObsDiagsQG : public util::Printable {
 public:
  typedef oops::Locations<QgObsTraits> Locations_;

  ObsDiagsQG(const ObsSpaceQG &, const Locations_ &, const oops::ObsVariables &) {}
  ~ObsDiagsQG() {}

// I/O
  void save(const std::string &) const {}

 private:
  void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------
}  // namespace qg

#endif  // QG_MODEL_OBSDIAGSQG_H_
