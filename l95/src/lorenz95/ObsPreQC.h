/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_OBSPREQC_H_
#define LORENZ95_OBSPREQC_H_

#include <memory>
#include <ostream>

#include "boost/shared_ptr.hpp"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace lorenz95 {
  class GomL95;
  class ObsTableView;
  template <typename DATATYPE> class ObsData1D;
  class ObsVec1D;

// Nothing to do here for the Lorenz model

class ObsPreQC : public util::Printable {
 public:
  ObsPreQC(ObsTableView &, const eckit::Configuration &,
           boost::shared_ptr<ObsData1D<int> >, boost::shared_ptr<ObsData1D<float> >);
  ~ObsPreQC() {}

  void priorFilter(const GomL95 &) const {}
  void postFilter(const ObsVec1D &) const {}

  const oops::Variables & requiredGeoVaLs() const {return novars_;}

 private:
  void print(std::ostream &) const {}
  const oops::Variables novars_;
};

}  // namespace lorenz95

#endif  // LORENZ95_OBSPREQC_H_
