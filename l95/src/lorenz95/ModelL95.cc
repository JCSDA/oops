/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "eckit/config/Configuration.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"

#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------
static oops::interface::ModelMaker<L95Traits, ModelL95> makermodel_("L95");

// -----------------------------------------------------------------------------

ModelL95::ModelL95(const Resolution & resol, const eckit::Configuration & config)
  : resol_(resol), f_(config.getDouble("f")),
    tstep_(util::Duration(config.getString("tstep"))),
    dt_(tstep_.toSeconds()/432000.0), vars_({"x"})
{
  oops::Log::info() << *this << std::endl;
  oops::Log::trace() << "ModelL95::ModelL95 created" << std::endl;
}

// -----------------------------------------------------------------------------
ModelL95::~ModelL95()
{
  oops::Log::trace() << "ModelL95::~ModelL95 destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelL95::initialize(StateL95 &) const {}
void ModelL95::finalize(StateL95 &) const {}
// -----------------------------------------------------------------------------
void ModelL95::step(StateL95 & xx, const ModelBias & bias) const {
  ModelTrajectory traj(false);
  this->stepRK(xx.getField(), bias, traj);
  xx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------

void ModelL95::stepRK(FieldL95 & xx, const ModelBias & bias,
                      ModelTrajectory & traj) const {
  FieldL95 dx(xx, false);
  FieldL95 zz(xx, false);
  FieldL95 dz(xx, false);

  zz = xx;
  traj.set(zz);
  this->tendencies(zz, bias.bias(), dz);
  dx = dz;

  zz = xx;
  zz.axpy(0.5, dz);
  traj.set(zz);
  this->tendencies(zz, bias.bias(), dz);
  dx.axpy(2.0, dz);

  zz = xx;
  zz.axpy(0.5, dz);
  traj.set(zz);
  this->tendencies(zz, bias.bias(), dz);
  dx.axpy(2.0, dz);

  zz = xx;
  zz += dz;
  traj.set(zz);
  this->tendencies(zz, bias.bias(), dz);
  dx += dz;

  const double zt = 1.0/6.0;
  xx.axpy(zt, dx);
}

// -----------------------------------------------------------------------------

#ifdef __INTEL_COMPILER
#pragma optimize("", off)
#endif
void ModelL95::tendencies(const FieldL95 & xx, const double & bias,
                          FieldL95 & dx) const {
  const int nn = resol_.npoints();
  // intel 19 is doing some agressive optimization of this loop that
  // is modifying the solution.
  for (int jj = 0; jj < nn; ++jj) {
    int jm2 = jj - 2;
    int jm1 = jj - 1;
    int jp1 = jj + 1;
    if (jm2 < 0) jm2 += nn;
    if (jm1 < 0) jm1 += nn;
    if (jp1 >= nn) jp1 -= nn;
    const double dxdt = -xx[jm2] * xx[jm1] + xx[jm1] * xx[jp1] - xx[jj] + f_ - bias;
    dx[jj] = dt_ * dxdt;
  }
}
#ifdef __INTEL_COMPILER
#pragma optimize("", on)
#endif

// -----------------------------------------------------------------------------

void ModelL95::print(std::ostream & os) const {
  os << "ModelL95: resol = " << resol_ << ", f = " << f_ << ", tstep = " << tstep_;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
