/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/ObsFilter.h"

#include <vector>

#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
ObsFilter::ObsFilter(const ObsTable & obsdb, const eckit::Configuration & conf,
           std::shared_ptr<ObsData1D<int> > qcflags, std::shared_ptr<ObsData1D<float> > obserr,
           const int iteration)
  : obsdb_(obsdb), qcflags_(qcflags), obserr_(obserr), novars_(),
    config_(conf.getSubConfiguration("obs filtering")),
    threshold_(config_.getFloat("threshold", 0.0)),
    inflation_(config_.getFloat("inflate obs error", -1.0)),
    bgCheck_(config_.has("threshold")),
    saveGeoVaLs_(config_.getBool("save geovals", false))
{
  oops::Log::info() << "config: " << conf << std::endl;
}

// -----------------------------------------------------------------------------
void ObsFilter::priorFilter(const GomL95 & gv) {
  if (saveGeoVaLs_) {
    gv.write(config_);
  }
}

// -----------------------------------------------------------------------------
void ObsFilter::postFilter(const GomL95 &, const ObsVec1D & hofx, const ObsVec1D &,
                                 const ObsDiags1D &) {
  if (bgCheck_) {
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
    oops::Log::info() << "ObsFilter::postFilter rejected = " << ireject
                      << ", inflated = " << inflate << std::endl;
  }
}

// -----------------------------------------------------------------------------
void ObsFilter::print(std::ostream & os) const {
  os << "L95 Background check with absolute threshold " << threshold_;
}

}  // namespace lorenz95

