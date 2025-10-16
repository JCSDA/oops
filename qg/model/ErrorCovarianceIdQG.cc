/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ErrorCovarianceIdQG.h"

#include <cmath>

#include "model/IncrementQG.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
ErrorCovarianceIdQG::ErrorCovarianceIdQG(const GeometryQG & resol, const oops::Variables & vars,
                                         const eckit::Configuration & conf,
                                         const StateQG &, const StateQG &)
{}
// -----------------------------------------------------------------------------
void ErrorCovarianceIdQG::multiply(const IncrementQG & dxin, IncrementQG & dxout) const {
  dxout = dxin;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceIdQG::inverseMultiply(const IncrementQG & dxin, IncrementQG & dxout) const {
  dxout = dxin;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceIdQG::randomize(IncrementQG & dx) const {
}
// -----------------------------------------------------------------------------
void ErrorCovarianceIdQG::print(std::ostream & os) const {
  os << "Identity error covariance";
}
// -----------------------------------------------------------------------------

}  // namespace qg
