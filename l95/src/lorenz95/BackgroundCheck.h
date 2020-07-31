/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_BACKGROUNDCHECK_H_
#define LORENZ95_BACKGROUNDCHECK_H_

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


/// Simple background check: all obs for which {|y-H(x)| < threshold} pass QC
class BackgroundCheck : public util::Printable {
 public:
  BackgroundCheck(const ObsTableView &, const eckit::Configuration &,
            boost::shared_ptr<ObsData1D<int> >, boost::shared_ptr<ObsData1D<float> >);

  void preProcess() const {}
  void priorFilter(const GomL95 &) const {}
  void postFilter(const ObsVec1D &, const ObsDiags1D &) const;

  oops::Variables requiredVars() const {return novars_;}
  oops::Variables requiredHdiagnostics() const {return novars_;}

 private:
  void print(std::ostream & os) const;

  const ObsTableView & obsdb_;
  const float threshold_;
  boost::shared_ptr<ObsData1D<int> > qcflags_;
  const oops::Variables novars_;
};

}  // namespace lorenz95

#endif  // LORENZ95_BACKGROUNDCHECK_H_
