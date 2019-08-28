/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_QCMANAGER_H_
#define LORENZ95_QCMANAGER_H_

#include <ostream>

#include "boost/shared_ptr.hpp"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace lorenz95 {
  class GomL95;
  template <typename DATATYPE> class ObsData1D;
  class ObsTableView;
  class ObsDiags1D;
  class ObsVec1D;

// Nothing to do here for the Lorenz model

class QCmanager : public util::Printable {
 public:
  QCmanager(const ObsTableView &, const eckit::Configuration &,
            boost::shared_ptr<ObsData1D<int> >, boost::shared_ptr<ObsData1D<float> >): novars_() {}
  ~QCmanager() {}

  void preProcess() const {}
  void priorFilter(const GomL95 &) const {}
  void postFilter(const ObsVec1D &, const ObsDiags1D &) const {}

  const oops::Variables & requiredGeoVaLs() const {return novars_;}
  const oops::Variables & requiredHdiagnostics() const {return novars_;}

 private:
  void print(std::ostream &) const {}
  const oops::Variables novars_;
};

}  // namespace lorenz95

#endif  // LORENZ95_QCMANAGER_H_
