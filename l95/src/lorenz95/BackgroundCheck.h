/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_BACKGROUNDCHECK_H_
#define LORENZ95_BACKGROUNDCHECK_H_

#include <memory>
#include <ostream>

#include "lorenz95/L95Traits.h"

#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class GomL95;
  template <typename DATATYPE> class ObsData1D;
  class ObsTable;
  class ObsDiags1D;
  class ObsVec1D;

/// Simple background check: all obs for which {|y-H(x)| < threshold} pass QC
class BackgroundCheck : public oops::interface::ObsFilterBase<L95ObsTraits> {
 public:
  BackgroundCheck(const ObsTable &, const eckit::Configuration &,
                  std::shared_ptr<ObsData1D<int> >, std::shared_ptr<ObsData1D<float> >);

  void preProcess() override {}
  void priorFilter(const GomL95 &) override {}
  void postFilter(const GomL95 &, const ObsVec1D &, const ObsVec1D &, const ObsDiags1D &) override;
  void checkFilterData(const oops::FilterStage filterStage) override {}

  oops::Variables requiredVars() const override {return novars_;}
  oops::Variables requiredHdiagnostics() const override {return novars_;}

 private:
  void print(std::ostream & os) const override;

  const ObsTable & obsdb_;
  std::shared_ptr<ObsData1D<int> > qcflags_;   // QC flags
  std::shared_ptr<ObsData1D<float> > obserr_;  // obs error stddev
  const oops::Variables novars_;
  const float threshold_;
  const float inflation_;
};

}  // namespace lorenz95

#endif  // LORENZ95_BACKGROUNDCHECK_H_
