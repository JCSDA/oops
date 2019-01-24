/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_OBSPREQC_H_
#define QG_MODEL_OBSPREQC_H_

#include <ostream>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace qg {
  class GomQG;
  class ObsSpaceQG;
  class ObsVecQG;

class ObsPreQC : public util::Printable {
 public:
  ObsPreQC(ObsSpaceQG &, const eckit::Configuration &);
  ~ObsPreQC() {}

  void priorFilter(const GomQG &) const {}
  void postFilter(const ObsVecQG &) const {}

  const oops::Variables & requiredGeoVaLs() const {return novars_;}

 private:
  void print(std::ostream &) const {}
  const oops::Variables novars_;
};

}  // namespace qg

#endif  // QG_MODEL_OBSPREQC_H_
