/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/FieldsQG.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/base/Variables.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "model/GomQG.h"
#include "model/LocationsQG.h"
#include "model/QgFortran.h"
#include "model/GeometryQG.h"
#include "model/VariablesQG.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const GeometryQG & geom, const oops::Variables & vars,
                   const util::DateTime & time):
  geom_(new GeometryQG(geom)), vars_(vars), time_(time)
{
  qg_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  qg_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  if (copy) {
    qg_field_copy_f90(keyFlds_, other.keyFlds_);
  } else {
    qg_field_zero_f90(keyFlds_);
  }
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  qg_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  qg_field_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const GeometryQG & geom)
  : geom_(new GeometryQG(geom)), vars_(other.vars_), time_(other.time_)
{
  qg_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  qg_field_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const oops::Variables & vars)
  : geom_(other.geom_), vars_(vars), time_(other.time_)
{
  qg_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  qg_field_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG::~FieldsQG() {
  qg_field_delete_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG & FieldsQG::operator=(const FieldsQG & rhs) {
  qg_field_copy_f90(keyFlds_, rhs.keyFlds_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
FieldsQG & FieldsQG::operator+=(const FieldsQG & rhs) {
  qg_field_self_add_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsQG & FieldsQG::operator-=(const FieldsQG & rhs) {
  qg_field_self_sub_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsQG & FieldsQG::operator*=(const double & zz) {
  qg_field_self_mul_f90(keyFlds_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void FieldsQG::zero() {
  qg_field_zero_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::zero(const util::DateTime & time) {
  qg_field_zero_f90(keyFlds_);
  time_ = time;
}
// -----------------------------------------------------------------------------
void FieldsQG::axpy(const double & zz, const FieldsQG & rhs) {
  qg_field_axpy_f90(keyFlds_, zz, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
double FieldsQG::dot_product_with(const FieldsQG & fld2) const {
  double zz;
  qg_field_dot_prod_f90(keyFlds_, fld2.keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsQG::schur_product_with(const FieldsQG & dx) {
    qg_field_self_schur_f90(keyFlds_, dx.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::random() {
  qg_field_random_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::dirac(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  qg_field_dirac_f90(keyFlds_, &conf);
}
// -----------------------------------------------------------------------------
void FieldsQG::interpolate(const LocationsQG & locs, const oops::Variables & vars,
                           GomQG & gom) const {
  const VariablesQG varqg(vars);
  qg_field_interp_tl_f90(keyFlds_, locs.toFortran(), varqg.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::interpolateTL(const LocationsQG & locs, const oops::Variables & vars,
                             GomQG & gom) const {
  const VariablesQG varqg(vars);
  qg_field_interp_tl_f90(keyFlds_, locs.toFortran(), varqg.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::interpolateAD(const LocationsQG & locs, const oops::Variables & vars,
                             const GomQG & gom) {
  const VariablesQG varqg(vars);
  qg_field_interp_ad_f90(keyFlds_, locs.toFortran(), varqg.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::changeResolution(const FieldsQG & other) {
  qg_field_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::add(const FieldsQG & rhs) {
  qg_field_add_incr_f90(keyFlds_, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::diff(const FieldsQG & x1, const FieldsQG & x2) {
  qg_field_diff_incr_f90(keyFlds_, x1.keyFlds_, x2.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::define(oops::UnstructuredGrid & ug) const {
  qg_field_define_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::convert_to(oops::UnstructuredGrid & ug) const {
  qg_field_convert_to_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::convert_from(const oops::UnstructuredGrid & ug) {
  qg_field_convert_from_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  qg_field_read_file_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
void FieldsQG::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  qg_field_write_file_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
double FieldsQG::norm() const {
  double zz = 0.0;
  qg_field_rms_f90(keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsQG::print(std::ostream & os) const {
  int nx = -1;
  int ny = -1;
  int nf = -1;
  int nb = -1;
  qg_field_sizes_f90(keyFlds_, nx, ny, nf, nb);
  os << std::endl << "  Resolution = " << nx << ", " << ny
     << ", Fields = " << nf << ", " << nb;
  nf += nb;
  std::vector<double> zstat(3*nf);
  qg_field_gpnorm_f90(keyFlds_, nf, zstat[0]);
  for (int jj = 0; jj < nf; ++jj) {
    os << std::endl << "  Min=" << zstat[3*jj]
       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2];
  }
}
// -----------------------------------------------------------------------------
bool FieldsQG::isForModel(bool nonlinear) const {
  int nx = -1;
  int ny = -1;
  int nf = -1;
  int nb = -1;
  qg_field_sizes_f90(keyFlds_, nx, ny, nf, nb);
  bool ok = (nf == 4);
  if (nonlinear) ok = ok && (nb == 2);
  return ok;
}
// -----------------------------------------------------------------------------
}  // namespace qg
