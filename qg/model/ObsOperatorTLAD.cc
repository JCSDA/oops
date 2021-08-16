/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsOperatorTLAD.h"

#include "eckit/config/Configuration.h"
#include "model/GomQG.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsOpBaseTLAD.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "oops/base/Variables.h"

namespace qg {

// -----------------------------------------------------------------------------

ObsOperatorTLAD::ObsOperatorTLAD(const ObsSpaceQG & os, const Parameters_ & params)
  : oper_(ObsOpTLADFactory::create(os, params.config))
{}

// -----------------------------------------------------------------------------

ObsOperatorTLAD::~ObsOperatorTLAD() {}

// -----------------------------------------------------------------------------

void ObsOperatorTLAD::setTrajectory(const GomQG & gvals, const ObsBias & bias) {
  oper_->setTrajectory(gvals, bias);
}

// -----------------------------------------------------------------------------

void ObsOperatorTLAD::simulateObsTL(const GomQG & gvals, ObsVecQG & yy,
                                    const ObsBiasIncrement & bias) const {
  oper_->simulateObsTL(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

void ObsOperatorTLAD::simulateObsAD(GomQG & gvals, const ObsVecQG & yy,
                                    ObsBiasIncrement & bias) const {
  oper_->simulateObsAD(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsOperatorTLAD::requiredVars() const {
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

void ObsOperatorTLAD::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace qg
