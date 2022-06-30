/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/BackgroundCheck.h"

#include <vector>

#include "lorenz95/L95Traits.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
static oops::interface::FilterMaker<L95ObsTraits, BackgroundCheck> makerBackgroundCheck_(
    "Background Check");

// -----------------------------------------------------------------------------
BackgroundCheck::BackgroundCheck(const ObsTable & obsdb, const Parameters_ & parameters,
           std::shared_ptr<ObsData1D<int> > qcflags, std::shared_ptr<ObsData1D<float> > obserr)
  : obsdb_(obsdb), options_(parameters), qcflags_(qcflags), obserr_(obserr), novars_()
{
}

// -----------------------------------------------------------------------------
void BackgroundCheck::postFilter(const GomL95 &,
                                 const ObsVec1D & hofx,
                                 const ObsVec1D &,
                                 const ObsDiags1D &) {
  std::vector<float> yobs;
  obsdb_.getdb("ObsValue", yobs);
  size_t inflate = 0;
  size_t ireject = 0;
  for (size_t jj = 0; jj < yobs.size(); ++jj) {
    if (std::abs(yobs[jj] - hofx[jj]) > options_.threshold) {
      // inflate obs error variance
      if (options_.inflation.value() != boost::none) {
        (*obserr_)[jj] *= *options_.inflation.value();
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
  os << "L95 Background check with absolute threshold " << options_.threshold;
}

}  // namespace lorenz95

