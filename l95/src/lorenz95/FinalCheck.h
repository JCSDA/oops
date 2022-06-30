/*
 * (C) Crown Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_FINALCHECK_H_
#define LORENZ95_FINALCHECK_H_

#include <memory>
#include <ostream>

#include "eckit/config/LocalConfiguration.h"

#include "lorenz95/L95Traits.h"

#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"

namespace lorenz95 {
  class GomL95;
  template <typename DATATYPE> class ObsData1D;
  class ObsTable;
  class ObsDiags1D;
  class ObsVec1D;

// Nothing to do here for the Lorenz model

class FinalCheck : public oops::interface::ObsFilterBase<L95ObsTraits> {
 public:
  FinalCheck(const ObsTable &, const eckit::Configuration &,
             std::shared_ptr<ObsData1D<int> >, std::shared_ptr<ObsData1D<float> >): novars_() {}

  void preProcess() override {}
  void priorFilter(const GomL95 &) override {}
  void postFilter(const GomL95 &,
                  const ObsVec1D &,
                  const ObsVec1D &,
                  const ObsDiags1D &) override {}
  void checkFilterData(const oops::FilterStage filterStage) override {}

  oops::Variables requiredVars() const override {return novars_;}
  oops::Variables requiredHdiagnostics() const override {return novars_;}

 private:
  void print(std::ostream &) const override {}
  const oops::Variables novars_;
};

}  // namespace lorenz95

#endif  // LORENZ95_FINALCHECK_H_
