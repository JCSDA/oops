/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/LocalizationMatrixQG.h"

#include "eckit/config/Configuration.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/QgFortran.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
LocalizationMatrixQG::LocalizationMatrixQG(const GeometryQG & resol,
                                           const oops::Variables & vars,
                                           const eckit::Configuration & config) {
  qg_error_covariance_setup_f90(keyLocal_, config, resol.toFortran());
}
// -----------------------------------------------------------------------------
LocalizationMatrixQG::~LocalizationMatrixQG() {
  qg_error_covariance_delete_f90(keyLocal_);
}
// -----------------------------------------------------------------------------
void LocalizationMatrixQG::randomize(IncrementQG & dx) const {
  qg_error_covariance_randomize_f90(keyLocal_, dx.fields().toFortran());
}
// -----------------------------------------------------------------------------
void LocalizationMatrixQG::multiply(IncrementQG & dx) const {
  IncrementQG dxtmp(dx);
  qg_error_covariance_mult_f90(keyLocal_, dxtmp.fields().toFortran(), dx.fields().toFortran());
}
// -----------------------------------------------------------------------------
void LocalizationMatrixQG::print(std::ostream & os) const {
  os << "LocalizationMatrixQG::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace qg
