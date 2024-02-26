/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/BackgroundCheck.h"

#include <vector>

#include "eckit/config/Configuration.h"

#include "lorenz95/L95Traits.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
static oops::interface::FilterMaker<L95ObsTraits, BackgroundCheck>
    makerBackgroundCheck_("Background Check");

// -----------------------------------------------------------------------------
BackgroundCheck::BackgroundCheck(const ObsTable & obsdb, const eckit::Configuration & conf,
           std::shared_ptr<ObsData1D<int> > qcflags, std::shared_ptr<ObsData1D<float> > obserr)
  : obsdb_(obsdb), qcflags_(qcflags), obserr_(obserr), novars_(),
    threshold_(conf.getFloat("threshold")), inflation_(conf.getFloat("inflate obs error", -1.0))
{}

// -----------------------------------------------------------------------------
void BackgroundCheck::postFilter(const GomL95 &, const ObsVec1D & hofx, const ObsVec1D &,
                                 const ObsDiags1D &) {
  std::vector<float> yobs;
  obsdb_.getdb("ObsValue", yobs);
  size_t inflate = 0;
  size_t ireject = 0;
  for (size_t jj = 0; jj < yobs.size(); ++jj) {
    if (std::abs(yobs[jj] - hofx[jj]) > threshold_) {
      // inflate obs error variance
      if (inflation_ > 0.0) {
        (*obserr_)[jj] *= inflation_;
        ++inflate;
      // or reject observation
      } else {
        (*qcflags_)[jj] = 1;
        ++ireject;
      }
    }
  }
  oops::Log::info() << "BackgroundCheck::postFilter rejected = " << ireject
                    << ", inflated = " << inflate << std::endl;
}

// -----------------------------------------------------------------------------
void BackgroundCheck::print(std::ostream & os) const {
  os << "L95 Background check with absolute threshold " << threshold_;
}

}  // namespace lorenz95

