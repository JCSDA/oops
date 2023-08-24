/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/TLML95.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/L95Traits.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"


namespace lorenz95 {
// -----------------------------------------------------------------------------
static oops::interface::LinearModelMaker<L95Traits, TLML95> makerTLML95_("L95TLM");
// -----------------------------------------------------------------------------
TLML95::TLML95(const Resolution & resol, const eckit::Configuration & tlConf)
  : resol_(resol), tstep_(util::Duration(tlConf.getString("tstep"))),
    dt_(tstep_.toSeconds()/432000.0), traj_(),
    lrmodel_(resol_, eckit::LocalConfiguration(tlConf, "trajectory")),
    vars_()
{
  oops::Log::info() << "TLML95: resol = " << resol_ << ", tstep = " << tstep_ << std::endl;
  oops::Log::trace() << "TLML95::TLML95 created" << std::endl;
}
// -----------------------------------------------------------------------------
TLML95::~TLML95() {
  for (trajIter jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
    delete jtra->second;
  }
  oops::Log::trace() << "TLML95::~TLML95 destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TLML95::setTrajectory(const StateL95 & xx, StateL95 &, const ModelBias & bias) {
  ASSERT(traj_.find(xx.validTime()) == traj_.end());
  ModelTrajectory * traj = new ModelTrajectory();
// Interpolate xx to xlr here
  FieldL95 zz(xx.getField());
  lrmodel_.stepRK(zz, bias, *traj);
  traj_[xx.validTime()] = traj;
}
// -----------------------------------------------------------------------------
const ModelTrajectory * TLML95::getTrajectory(const util::DateTime & tt) const {
  trajICst itra = traj_.find(tt);
  if (itra == traj_.end()) {
    oops::Log::error() << "TLML95: trajectory not available at time " << tt << std::endl;
    ABORT("TLML95: trajectory not available");
  }
  return itra->second;
}
// -----------------------------------------------------------------------------
/// Run TLM and its adjoint
// -----------------------------------------------------------------------------
void TLML95::initializeTL(IncrementL95 &) const {}
void TLML95::finalizeTL(IncrementL95 &) const {}
void TLML95::initializeAD(IncrementL95 &) const {}
void TLML95::finalizeAD(IncrementL95 &) const {}
// -----------------------------------------------------------------------------
void TLML95::stepTL(IncrementL95 & xx, const ModelBiasCorrection & bias) const {
  FieldL95 dx(xx.getField(), false);
  FieldL95 zz(xx.getField(), false);
  FieldL95 dz(xx.getField(), false);
  const ModelTrajectory * traj = this->getTrajectory(xx.validTime());

  zz = xx.getField();
  this->tendenciesTL(zz, bias.bias(), traj->get(1), dz);
  dx = dz;

  zz = xx.getField();
  zz.axpy(0.5, dz);
  this->tendenciesTL(zz, bias.bias(), traj->get(2), dz);
  dx.axpy(2.0, dz);

  zz = xx.getField();
  zz.axpy(0.5, dz);
  this->tendenciesTL(zz, bias.bias(), traj->get(3), dz);
  dx.axpy(2.0, dz);

  zz = xx.getField();
  zz += dz;
  this->tendenciesTL(zz, bias.bias(), traj->get(4), dz);
  dx += dz;

  const double zt = 1.0/6.0;
  xx.getField().axpy(zt, dx);
  xx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------
void TLML95::stepAD(IncrementL95 & xx, ModelBiasCorrection & bias) const {
  FieldL95 dx(xx.getField(), false);
  FieldL95 zz(xx.getField(), false);
  FieldL95 dz(xx.getField(), false);

  xx.validTime() -= tstep_;
  const ModelTrajectory * traj = this->getTrajectory(xx.validTime());

  const double zt = 1.0/6.0;
  dx = xx.getField();
  dx *= zt;

  dz = dx;
  this->tendenciesAD(zz, bias.bias(), traj->get(4), dz);
  xx.getField() += zz;
  dz = zz;

  dz.axpy(2.0, dx);
  this->tendenciesAD(zz, bias.bias(), traj->get(3), dz);
  xx.getField() += zz;
  dz = zz;
  dz *= 0.5;

  dz.axpy(2.0, dx);
  this->tendenciesAD(zz, bias.bias(), traj->get(2), dz);
  xx.getField() += zz;
  dz = zz;
  dz *= 0.5;

  dz += dx;
  this->tendenciesAD(zz, bias.bias(), traj->get(1), dz);
  xx.getField() += zz;
}
// -----------------------------------------------------------------------------
// intel 19 tries to aggressive optimize these functions in a way that leads
// to memory violations.  So, turn off optimizations.
#ifdef __INTEL_COMPILER
#pragma optimize("", off)
#endif
void TLML95::tendenciesTL(const FieldL95 & xx, const double & bias,
                          const FieldL95 & xtraj, FieldL95 & dx) const {
  const int nn = resol_.npoints();
  for (int jj = 0; jj < nn; ++jj) {
    int jm2 = jj - 2;
    int jm1 = jj - 1;
    int jp1 = jj + 1;
    if (jm2 < 0) jm2 += nn;
    if (jm1 < 0) jm1 += nn;
    if (jp1 >= nn) jp1 -= nn;
    const double dxdt = - xx[jm2] * xtraj[jm1] - xtraj[jm2] * xx[jm1]
                        + xx[jm1] * xtraj[jp1] + xtraj[jm1] * xx[jp1]
                        - xx[jj] - bias;
    dx[jj] = dt_ * dxdt;
  }
}
// -----------------------------------------------------------------------------
void TLML95::tendenciesAD(FieldL95 & xx, double & bias,
                          const FieldL95 & xtraj, const FieldL95 & dx) const {
  const int nn = resol_.npoints();
  xx.zero();
  for (int jj = 0; jj < nn; ++jj) {
    int jm2 = jj - 2;
    int jm1 = jj - 1;
    int jp1 = jj + 1;
    if (jm2 < 0) jm2 += nn;
    if (jm1 < 0) jm1 += nn;
    if (jp1 >= nn) jp1 -= nn;
    const double dxdt = dt_ * dx[jj];
    xx[jm2] -= dxdt * xtraj[jm1];
    xx[jm1] -= dxdt * xtraj[jm2];
    xx[jm1] += dxdt * xtraj[jp1];
    xx[jp1] += dxdt * xtraj[jm1];
    xx[jj] -= dxdt;
    bias -= dxdt;
  }
}
#ifdef __INTEL_COMPILER
#pragma optimize("", on)
#endif
// -----------------------------------------------------------------------------
void TLML95::print(std::ostream & os) const {
  os << "TLML95: resol = " << resol_ << ", tstep = " << tstep_ << std::endl;
  os << "L95 Model Trajectory, nstep=" << traj_.size() << std::endl;
  typedef std::map< util::DateTime, ModelTrajectory * >::const_iterator trajICst;
  if (traj_.size() > 0) {
    os << "L95 Model Trajectory: times are:";
    for (trajICst jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
      os << "  " << jtra->first;
    }
  }
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
