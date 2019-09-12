/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_QCMANAGER_H_
#define QG_MODEL_QCMANAGER_H_

#include <ostream>

#include "boost/shared_ptr.hpp"

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace qg {
  class GomQG;
  template <typename DATATYPE> class ObsDataQG;
  class ObsDiagsQG;
  class ObsSpaceQG;
  class ObsVecQG;

class QCmanager : public util::Printable {
 public:
  QCmanager(const ObsSpaceQG &, const eckit::Configuration &,
            boost::shared_ptr<ObsDataQG<int> >, boost::shared_ptr<ObsDataQG<float> >): novars_() {}
  ~QCmanager() {}

  void preProcess() const {}
  void priorFilter(const GomQG &) const {}
  void postFilter(const ObsVecQG &, const ObsDiagsQG &) const {}

  oops::Variables requiredGeoVaLs() const {return novars_;}
  oops::Variables requiredHdiagnostics() const {return novars_;}

 private:
  void print(std::ostream &) const {}
  const oops::Variables novars_;
};

}  // namespace qg

#endif  // QG_MODEL_QCMANAGER_H_
