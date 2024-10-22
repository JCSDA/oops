/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSFILTER_H_
#define LORENZ95_OBSFILTER_H_

#include <memory>
#include <ostream>

#include "lorenz95/L95Traits.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class GomL95;
  template <typename DATATYPE> class ObsData1D;
  class ObsTable;
  class ObsDiags1D;
  class ObsVec1D;

/// ObsFilter includes optional:
/// - simple background check: all obs for which {|y-H(x)| < threshold} pass QC
/// - geovals saver
class ObsFilter : public util::Printable {
 public:
  ObsFilter(const ObsTable &, const eckit::Configuration &,
                  std::shared_ptr<ObsData1D<int> >, std::shared_ptr<ObsData1D<float> >,
                  const int iteration = 0);

  void preProcess() {}
  void priorFilter(const GomL95 &);
  void postFilter(const GomL95 &, const ObsVec1D &, const ObsVec1D &, const ObsDiags1D &);

  oops::Variables requiredVars() const {return novars_;}
  oops::ObsVariables requiredHdiagnostics() const {return noobsvars_;}

 private:
  void print(std::ostream & os) const override;

  const ObsTable & obsdb_;
  std::shared_ptr<ObsData1D<int> > qcflags_;   // QC flags
  std::shared_ptr<ObsData1D<float> > obserr_;  // obs error stddev
  const oops::Variables novars_;
  const oops::ObsVariables noobsvars_;
  const eckit::LocalConfiguration config_;
  const float threshold_;
  const float inflation_;
  const bool bgCheck_;
  const bool saveGeoVaLs_;
};

}  // namespace lorenz95

#endif  // LORENZ95_OBSFILTER_H_
