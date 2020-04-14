/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ErrorStdDevQG.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/StateQG.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
ErrorStdDevQG::ErrorStdDevQG(const StateQG &, const StateQG &,
                             const GeometryQG &, const eckit::Configuration & conf)
  : keyFtnConfig_(0)
{
  qg_error_stddev_setup_f90(keyFtnConfig_, conf);
  oops::Log::trace() << "ErrorStdDevQG created" << std::endl;
}
// -----------------------------------------------------------------------------
ErrorStdDevQG::~ErrorStdDevQG() {
  qg_error_stddev_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ErrorStdDevQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ErrorStdDevQG::multiply(const IncrementQG & dxin, IncrementQG & dxout) const {
  qg_error_stddev_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
                           dxout.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ErrorStdDevQG::multiplyInverse(const IncrementQG & dxin, IncrementQG & dxout) const {
  qg_error_stddev_inv_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
                              dxout.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ErrorStdDevQG::multiplyAD(const IncrementQG & dxin, IncrementQG & dxout) const {
  qg_error_stddev_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
                           dxout.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ErrorStdDevQG::multiplyInverseAD(const IncrementQG & dxin, IncrementQG & dxout) const {
  qg_error_stddev_inv_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
                              dxout.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ErrorStdDevQG::print(std::ostream & os) const {
  os << "QG Background Error Standard Deviations";
}
// -----------------------------------------------------------------------------
}  // namespace qg

