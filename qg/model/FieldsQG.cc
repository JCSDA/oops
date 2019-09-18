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
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "model/GeometryQG.h"
#include "model/GomQG.h"
#include "model/LocationsQG.h"
#include "model/QgFortran.h"
#include "model/VariablesQG.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const GeometryQG & geom, const oops::Variables & vars,
                   const bool & lbc, const util::DateTime & time):
  geom_(new GeometryQG(geom)), vars_(vars), lbc_(lbc), time_(time)
{
  qg_fields_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran(), lbc_);
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), lbc_(other.lbc_), time_(other.time_)
{
  qg_fields_create_from_other_f90(keyFlds_, other.keyFlds_);
  if (copy) {
    qg_fields_copy_f90(keyFlds_, other.keyFlds_);
  }
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other)
  : geom_(other.geom_), vars_(other.vars_), lbc_(other.lbc_), time_(other.time_)
{
  qg_fields_create_from_other_f90(keyFlds_, other.keyFlds_);
  qg_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const GeometryQG & geom)
  : geom_(new GeometryQG(geom)), vars_(other.vars_), lbc_(other.lbc_), time_(other.time_)
{
  qg_fields_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran(), lbc_);
  qg_fields_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const oops::Variables & vars)
  : geom_(other.geom_), vars_(vars), lbc_(other.lbc_), time_(other.time_)
{
// TODO(Benjamin): delete that ?
  qg_fields_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran(), lbc_);
  qg_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG::~FieldsQG() {
  qg_fields_delete_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG & FieldsQG::operator=(const FieldsQG & rhs) {
  qg_fields_copy_f90(keyFlds_, rhs.keyFlds_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
FieldsQG & FieldsQG::operator+=(const FieldsQG & rhs) {
  qg_fields_self_add_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsQG & FieldsQG::operator-=(const FieldsQG & rhs) {
  qg_fields_self_sub_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsQG & FieldsQG::operator*=(const double & zz) {
  qg_fields_self_mul_f90(keyFlds_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void FieldsQG::zero() {
  qg_fields_zero_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::zero(const util::DateTime & time) {
  qg_fields_zero_f90(keyFlds_);
  time_ = time;
}
// -----------------------------------------------------------------------------
void FieldsQG::axpy(const double & zz, const FieldsQG & rhs) {
  qg_fields_axpy_f90(keyFlds_, zz, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
double FieldsQG::dot_product_with(const FieldsQG & fld2) const {
  double zz;
  qg_fields_dot_prod_f90(keyFlds_, fld2.keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsQG::schur_product_with(const FieldsQG & dx) {
    qg_fields_self_schur_f90(keyFlds_, dx.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::random() {
  qg_fields_random_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::dirac(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  qg_fields_dirac_f90(keyFlds_, &conf);
}
// -----------------------------------------------------------------------------
void FieldsQG::getValues(const LocationsQG & locs, const oops::Variables & vars,
                         GomQG & gom) const {
  const VariablesQG varqg(vars);
  qg_fields_interp_f90(keyFlds_, locs.toFortran(), varqg.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::getValuesTL(const LocationsQG & locs, const oops::Variables & vars,
                           GomQG & gom) const {
  const VariablesQG varqg(vars);
  qg_fields_interp_tl_f90(keyFlds_, locs.toFortran(), varqg.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::getValuesAD(const LocationsQG & locs, const oops::Variables & vars,
                           const GomQG & gom) {
  const VariablesQG varqg(vars);
  qg_fields_interp_ad_f90(keyFlds_, locs.toFortran(), varqg.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::changeResolution(const FieldsQG & other) {
  qg_fields_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::add(const FieldsQG & rhs) {
  FieldsQG rhs_myres(rhs, *geom_);
  qg_fields_add_incr_f90(keyFlds_, rhs_myres.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::diff(const FieldsQG & x1, const FieldsQG & x2) {
  FieldsQG x1_myres(x1, *geom_);
  FieldsQG x2_myres(x2, *geom_);
  qg_fields_diff_incr_f90(keyFlds_, x1_myres.keyFlds_, x2_myres.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsQG::ug_coord(oops::UnstructuredGrid & ug) const {
  qg_fields_ug_coord_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsQG::field_to_ug(oops::UnstructuredGrid & ug, const int & its) const {
  qg_fields_field_to_ug_f90(keyFlds_, ug.toFortran(), its);
}
// -----------------------------------------------------------------------------
void FieldsQG::field_from_ug(const oops::UnstructuredGrid & ug, const int & its) {
  qg_fields_field_from_ug_f90(keyFlds_, ug.toFortran(), its);
}
// -----------------------------------------------------------------------------
void FieldsQG::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  qg_fields_read_file_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
void FieldsQG::analytic_init(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  qg_fields_analytic_init_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
void FieldsQG::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  qg_fields_write_file_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
double FieldsQG::norm() const {
  double zz = 0.0;
  qg_fields_rms_f90(keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsQG::print(std::ostream & os) const {
  int nx, ny, nz, nb, lq, lbc;
  qg_fields_sizes_f90(keyFlds_, nx, ny, nz, nb);
  qg_fields_vars_f90(keyFlds_, lq, lbc);
  os << std::endl << "  Resolution = " << nx << ", " << ny << ", " << nz;
  if (lq == 1) {
    os << std::endl << "  Variable = potential vorticity";
  } else {
    os << std::endl << "  Variable = streamfunction";
  }
  if (lbc == 1) {
    os << std::endl << "  Boundary conditions are activated";
  } else {
    os << std::endl << "  Boundary conditions are not activated";
  }
  std::vector<double> zstat(4*(1+nb));
  qg_fields_gpnorm_f90(keyFlds_, nb, zstat[0]);
  for (int jj = 0; jj < 1+nb; ++jj) {
    std::ios_base::fmtflags f(os.flags());
    os << std::endl << "  Scaling=" << std::setprecision(4) << std::setw(7) << zstat[4*jj]
       << ", Min=" << std::fixed << std::setprecision(4) << std::setw(12) << zstat[4*jj+1]
       << ", Max=" << std::fixed << std::setprecision(4) << std::setw(12) <<zstat[4*jj+2]
       << ", RMS=" << std::fixed << std::setprecision(4) << std::setw(12) <<zstat[4*jj+3];
    os.flags(f);
  }
}
// -----------------------------------------------------------------------------
bool FieldsQG::isForModel(const bool & nonlinear) const {
  int nx, ny, nz, nb;
  qg_fields_sizes_f90(keyFlds_, nx, ny, nz, nb);
  bool ok = true;
  if (nonlinear) ok = (nb == 2);
  return ok;
}
// -----------------------------------------------------------------------------
oops::GridPoint FieldsQG::getPoint(const GeometryQGIterator & iter) const {
  int nx, ny, nz, nb;
  qg_fields_sizes_f90(keyFlds_, nx, ny, nz, nb);
  std::vector<int> varlens(1, nz);
  std::vector<double> values(nz);
  qg_fields_getpoint_f90(keyFlds_, iter.toFortran(), nz, values[0]);
  return oops::GridPoint(vars_.toOopsVariables(), values, varlens);
}
// -----------------------------------------------------------------------------
void FieldsQG::setPoint(const oops::GridPoint & x, const GeometryQGIterator & iter) {
  const std::vector<double> vals = x.getVals();
  qg_fields_setpoint_f90(keyFlds_, iter.toFortran(), vals.size(), vals[0]);
}
// -----------------------------------------------------------------------------
size_t FieldsQG::serialSize() const {
  size_t nn = 0;
  int nx, ny, nz, nb;
  qg_fields_sizes_f90(keyFlds_, nx, ny, nz, nb);
  nn += nx * ny * nz + nb * nx * nz;
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void FieldsQG::serialize(std::vector<double> & vect)  const {
  int size_fld = this->serialSize() - 2;

  // Allocate space for fld, xb and qb
  std::vector<double> v_fld(size_fld, 0);

  // Serialize the field
  qg_fields_serialize_f90(keyFlds_, size_fld, v_fld.data());
  vect.insert(vect.end(), v_fld.begin(), v_fld.end());

  // Serialize the date and time
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void FieldsQG::deserialize(const std::vector<double> & vect, size_t & index) {
  qg_fields_deserialize_f90(keyFlds_, vect.size(), vect.data(), index);
  time_.deserialize(vect, index);
}
// -----------------------------------------------------------------------------
}  // namespace qg
