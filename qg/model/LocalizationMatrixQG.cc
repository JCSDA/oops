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

#include "model/QgFortran.h"
#include "model/IncrementQG.h"
#include "model/StateQG.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
LocalizationMatrixQG::LocalizationMatrixQG(const StateQG & xb,
                                           const eckit::Configuration & config) {
  const eckit::Configuration * configc = &config;
  qg_localization_setup_f90(keyFtnConfig_, &configc, xb.geometry()->toFortran());
}
// -----------------------------------------------------------------------------
LocalizationMatrixQG::~LocalizationMatrixQG() {
  qg_localization_delete_f90(keyFtnConfig_);
}
// -----------------------------------------------------------------------------
void LocalizationMatrixQG::multiply(IncrementQG & dx) const {
  qg_localization_mult_f90(keyFtnConfig_, dx.fields().toFortran());
}
// -----------------------------------------------------------------------------
void LocalizationMatrixQG::print(std::ostream & os) const {
  os << "LocalizationMatrixQG::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace qg
