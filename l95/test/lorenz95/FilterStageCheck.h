/*
 * (C) Copyright 2022 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef L95_TEST_LORENZ95_FILTERSTAGECHECK_H_
#define L95_TEST_LORENZ95_FILTERSTAGECHECK_H_

#include <memory>
#include <ostream>

#include "lorenz95/L95Traits.h"

#include "oops/base/Variables.h"
#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace test {

/// Parameters for L95 FilterStageCheck
class FilterStageCheckParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(FilterStageCheckParameters, ObsFilterParametersBase)
};

/// Filter stage check: reject different obs depending on whether the filter
/// is run at the pre, prior or post stage.
class FilterStageCheck : public oops::interface::ObsFilterBase<lorenz95::L95ObsTraits> {
 public:
  typedef FilterStageCheckParameters Parameters_;

  FilterStageCheck(const lorenz95::ObsTable &, const Parameters_ &,
                   std::shared_ptr<lorenz95::ObsData1D<int> >,
                   std::shared_ptr<lorenz95::ObsData1D<float> >);

  void preProcess() override;
  void priorFilter(const lorenz95::GomL95 &) override;
  void postFilter(const lorenz95::GomL95 &,
                  const lorenz95::ObsVec1D &,
                  const lorenz95::ObsVec1D &,
                  const lorenz95::ObsDiags1D &) override;
  void checkFilterData(const oops::FilterStage filterStage) override {}

  oops::Variables requiredVars() const override {return novars_;}
  oops::Variables requiredHdiagnostics() const override {return novars_;}

 private:
  void print(std::ostream & os) const override;

  const lorenz95::ObsTable & obsdb_;
  Parameters_ options_;
  std::shared_ptr<lorenz95::ObsData1D<int> > qcflags_;   // QC flags
  std::shared_ptr<lorenz95::ObsData1D<float> > obserr_;  // obs error stddev
  const oops::Variables novars_;
};

}  // namespace test

#endif  // L95_TEST_LORENZ95_FILTERSTAGECHECK_H_
