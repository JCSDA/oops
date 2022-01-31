/*
 * (C) Copyright 2022 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "../../../l95/test/lorenz95/FilterStageCheck.h"

#include <vector>

#include "lorenz95/L95Traits.h"

// -----------------------------------------------------------------------------
namespace test {
// -----------------------------------------------------------------------------
static oops::interface::FilterMaker<lorenz95::L95ObsTraits, FilterStageCheck>
makerFilterStageCheck_("Filter Stage Check");

// -----------------------------------------------------------------------------
FilterStageCheck::FilterStageCheck(const lorenz95::ObsTable & obsdb, const Parameters_ & parameters,
                                   std::shared_ptr<lorenz95::ObsData1D<int> > qcflags,
                                   std::shared_ptr<lorenz95::ObsData1D<float> > obserr)
  : obsdb_(obsdb), options_(parameters), qcflags_(qcflags), obserr_(obserr), novars_()
{
}

// -----------------------------------------------------------------------------
void FilterStageCheck::preProcess() {
  // Keep all observations
  oops::Log::info() << "FilterStageCheck::preProcess rejected nothing" << std::endl;
}

// -----------------------------------------------------------------------------
void FilterStageCheck::priorFilter(const lorenz95::GomL95 &) {
  std::vector<float> yobs;
  obsdb_.getdb("ObsValue", yobs);
  size_t ireject = 0;
  // Reject one third of the observations
  for (size_t jj = 0; jj < yobs.size() / 3; ++jj) {
    (*qcflags_)[jj] = 1;
    ++ireject;
  }
  oops::Log::info() << "FilterStageCheck::postFilter rejected = " << ireject << std::endl;
}

// -----------------------------------------------------------------------------
void FilterStageCheck::postFilter(const lorenz95::GomL95 &,
                                  const lorenz95::ObsVec1D &,
                                  const lorenz95::ObsVec1D &,
                                  const lorenz95::ObsDiags1D &) {
  std::vector<float> yobs;
  obsdb_.getdb("ObsValue", yobs);
  size_t ireject = 0;
  // Reject half of the observations
  for (size_t jj = 0; jj < yobs.size() / 2; ++jj) {
    (*qcflags_)[jj] = 1;
    ++ireject;
  }
  oops::Log::info() << "FilterStageCheck::postFilter rejected = " << ireject << std::endl;
}

// -----------------------------------------------------------------------------
void FilterStageCheck::print(std::ostream & os) const {
  os << "L95 Filter Stage Check" << std::endl;
}

}  // namespace test

