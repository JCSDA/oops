/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsOperatorQG.h"

#include "eckit/config/Configuration.h"
#include "model/GomQG.h"
#include "model/LocationsQG.h"
#include "model/ObsBias.h"
#include "model/ObsOpBaseQG.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

namespace qg {

// -----------------------------------------------------------------------------

ObsOperatorQG::ObsOperatorQG(const ObsSpaceQG & os, const eckit::Configuration & conf)
  : oper_(ObsOpFactory::create(os, conf)), obsv_({"zz"}, "")
{}

// -----------------------------------------------------------------------------

ObsOperatorQG::~ObsOperatorQG() {}

// -----------------------------------------------------------------------------

void ObsOperatorQG::simulateObs(const GomQG & gvals, ObsVecQG & yy, const ObsBias & bias) const {
  oper_->simulateObs(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsOperatorQG::variables() const {
  return oper_->variables();
}

// -----------------------------------------------------------------------------

LocationsQG * ObsOperatorQG::locations(const util::DateTime & t1, const util::DateTime & t2) const {
  return oper_->locations(t1, t2);
}

// -----------------------------------------------------------------------------

const std::string & ObsOperatorQG::obstype() const {
  return oper_->obstype();
}

// -----------------------------------------------------------------------------

void ObsOperatorQG::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace qg
