/*
 * (C) Crown Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_FINALCHECK_H_
#define QG_MODEL_FINALCHECK_H_

#include <memory>
#include <ostream>

#include "eckit/config/LocalConfiguration.h"

#include "model/QgTraits.h"

#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/Printable.h"

namespace qg {
  class GomQG;
  template <typename DATATYPE> class ObsDataQG;
  class ObsDiagsQG;
  class ObsSpaceQG;
  class ObsVecQG;

class FinalCheck : public oops::interface::ObsFilterBase<QgObsTraits> {
 public:
  FinalCheck(const ObsSpaceQG &, const eckit::Configuration &,
             std::shared_ptr<ObsDataQG<int> >, std::shared_ptr<ObsDataQG<float> >): novars_() {}

  void preProcess() override {}
  void priorFilter(const GomQG &) override {}
  void postFilter(const GomQG &,
                  const ObsVecQG &,
                  const ObsVecQG &,
                  const ObsDiagsQG &) override {}
  void checkFilterData(const oops::FilterStage filterStage) override {}

  oops::Variables requiredVars() const override {return novars_;}
  oops::Variables requiredHdiagnostics() const override {return novars_;}

 private:
  void print(std::ostream &) const override {}
  const oops::Variables novars_;
};

}  // namespace qg

#endif  // QG_MODEL_FINALCHECK_H_
