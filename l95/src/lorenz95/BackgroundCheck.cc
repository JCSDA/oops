/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/BackgroundCheck.h"

#include <vector>

#include "lorenz95/L95Traits.h"
#include "oops/interface/ObsFilter.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
static oops::FilterMaker<L95ObsTraits,
       oops::ObsFilter<L95ObsTraits, BackgroundCheck> > makerBackgroundCheck_("Background Check");

// -----------------------------------------------------------------------------
BackgroundCheck::BackgroundCheck(const ObsTableView & obsdb, const eckit::Configuration & conf,
        boost::shared_ptr<ObsData1D<int> > qcflags, boost::shared_ptr<ObsData1D<float> >)
  : obsdb_(obsdb), threshold_(conf.getFloat("threshold")), qcflags_(qcflags), novars_()
{
}

// -----------------------------------------------------------------------------
void BackgroundCheck::postFilter(const ObsVec1D & hofx, const ObsDiags1D &) const {
  std::vector<float> yobs;
  obsdb_.getdb("ObsValue", yobs);
  for (size_t jj = 0; jj < yobs.size(); ++jj) {
    if (std::abs(yobs[jj] - hofx[jj]) > threshold_) {
      (*qcflags_)[jj] = 1;
    }
  }
}

// -----------------------------------------------------------------------------
void BackgroundCheck::print(std::ostream & os) const {
  os << "L95 Background check with absolute threshold " << threshold_ << std::endl;
}

}  // namespace lorenz95

