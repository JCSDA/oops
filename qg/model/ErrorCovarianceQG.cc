/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ErrorCovarianceQG.h"

#include <cmath>

#include "eckit/config/Configuration.h"
#include "model/FieldsQG.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/QgFortran.h"
#include "model/StateQG.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------

ErrorCovarianceQG::ErrorCovarianceQG(const GeometryQG & resol, const oops::Variables &,
                                     const eckit::Configuration & conf,
                                     const StateQG &, const StateQG &) {
  time_ = util::DateTime(conf.getString("date"));
  const eckit::Configuration * configc = &conf;
  qg_b_setup_f90(keyFtnConfig_, &configc, resol.toFortran());
  oops::Log::trace() << "ErrorCovarianceQG created" << std::endl;
}

// -----------------------------------------------------------------------------

ErrorCovarianceQG::~ErrorCovarianceQG() {
  qg_b_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ErrorCovarianceQG destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ErrorCovarianceQG::linearize(const StateQG &, const GeometryQG & resol) {
  geom_.reset(new GeometryQG(resol));
}

// -----------------------------------------------------------------------------

void ErrorCovarianceQG::multiply(const IncrementQG & dxin, IncrementQG & dxout) const {
  qg_b_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
                               dxout.fields().toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceQG::inverseMultiply(const IncrementQG & dxin, IncrementQG & dxout) const {
  qg_b_invmult_f90(keyFtnConfig_, dxin.fields().toFortran(),
                               dxout.fields().toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceQG::randomize(IncrementQG & dx) const {
  qg_b_randomize_f90(keyFtnConfig_, dx.fields().toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceQG::print(std::ostream & os) const {
  os << "ErrorCovarianceQG::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace qg
