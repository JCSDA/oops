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

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "model/GeometryQG.h"
#include "model/GomQG.h"
#include "model/QgFortran.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const GeometryQG & geom, const oops::Variables & vars,
                   const bool & lbc, const util::DateTime & time):
  geom_(new GeometryQG(geom)), vars_(vars), lbc_(lbc), time_(time)
{
  qg_fields_create_f90(keyFlds_, geom_->toFortran(), vars_, lbc_);
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), lbc_(other.lbc_), time_(other.time_)
{
  qg_fields_create_from_other_f90(keyFlds_, other.keyFlds_, geom_->toFortran());
  if (copy) {
    qg_fields_copy_f90(keyFlds_, other.keyFlds_);
  }
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other)
  : geom_(other.geom_), vars_(other.vars_), lbc_(other.lbc_), time_(other.time_)
{
  qg_fields_create_from_other_f90(keyFlds_, other.keyFlds_, geom_->toFortran());
  qg_fields_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const GeometryQG & geom, const bool ad)
  : geom_(new GeometryQG(geom)), vars_(other.vars_), lbc_(other.lbc_), time_(other.time_)
{
  qg_fields_create_f90(keyFlds_, geom_->toFortran(), vars_, lbc_);
  if (ad) {
    qg_fields_change_resol_ad_f90(keyFlds_, other.keyFlds_);
  } else {
    qg_fields_change_resol_f90(keyFlds_, other.keyFlds_);
  }
}
// -----------------------------------------------------------------------------
FieldsQG::FieldsQG(const FieldsQG & other, const oops::Variables & vars)
  : geom_(other.geom_), vars_(vars), lbc_(other.lbc_), time_(other.time_)
{
// TODO(Benjamin): delete that ?
  qg_fields_create_f90(keyFlds_, geom_->toFortran(), vars_, lbc_);
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
void FieldsQG::ones() {
  qg_fields_ones_f90(keyFlds_);
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
  qg_fields_dirac_f90(keyFlds_, config);
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
void FieldsQG::toFieldSet(atlas::FieldSet & afieldset) const {
  qg_fields_to_fieldset_f90(keyFlds_, afieldset.get());
}
// -----------------------------------------------------------------------------
void FieldsQG::fromFieldSet(const atlas::FieldSet & afieldset) {
  qg_fields_from_fieldset_f90(keyFlds_, afieldset.get());
}
// -----------------------------------------------------------------------------
void FieldsQG::read(const eckit::Configuration & config) {
  qg_fields_read_file_f90(keyFlds_, config, time_);
}
// -----------------------------------------------------------------------------
void FieldsQG::analytic_init(const eckit::Configuration & config) {
  qg_fields_analytic_init_f90(keyFlds_, config, time_);
}
// -----------------------------------------------------------------------------
void FieldsQG::write(const eckit::Configuration & config) const {
  qg_fields_write_file_f90(keyFlds_, config, time_);
}
// -----------------------------------------------------------------------------
double FieldsQG::norm() const {
  double zz = 0.0;
  qg_fields_rms_f90(keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsQG::print(std::ostream & os) const {
  // Resolution
  int nx, ny, nz;
  qg_fields_sizes_f90(keyFlds_, nx, ny, nz);
  os << std::endl << "  Resolution = " << nx << ", " << ny << ", " << nz;

  // Min, max, RMS of fields
  std::vector<std::string> var{"Streamfunction         :",
                               "Potential vorticity    :",
                               "U-wind                 :",
                               "V-wind                 :",
                               "Streamfunction LBC     :",
                               "Potential vorticity LBC:"};
  std::vector<int> vpresent(6);
  std::vector<double> vmin(6);
  std::vector<double> vmax(6);
  std::vector<double> vrms(6);
  qg_fields_gpnorm_f90(keyFlds_, vpresent.data(), vmin.data(), vmax.data(), vrms.data());
  for (int jj = 0; jj < 6; ++jj) {
    if (vpresent[jj] == 1) {
      std::ios_base::fmtflags f(os.flags());
      os << std::endl << "  " << var[jj]
       << "  Min=" << vmin[jj]
       << ", Max=" << vmax[jj]
       << ", RMS=" << vrms[jj];
      os.flags(f);
    }
  }
}
// -----------------------------------------------------------------------------
bool FieldsQG::isForModel(const bool & nonlinear) const {
  bool ok = true;
  if (nonlinear) {
    int lbc;
    qg_fields_lbc_f90(keyFlds_, lbc);
    ok = (lbc == 1);
  }
  return ok;
}
// -----------------------------------------------------------------------------
oops::LocalIncrement FieldsQG::getLocal(const GeometryQGIterator & iter) const {
  int nx, ny, nz;
  qg_fields_sizes_f90(keyFlds_, nx, ny, nz);
  std::vector<int> varlens(vars_.size());
  for (unsigned int ii = 0; ii < vars_.size(); ii++) {
    varlens[ii] = nz;
  }
  int lenvalues = std::accumulate(varlens.begin(), varlens.end(), 0);
  std::vector<double> values(lenvalues);
  qg_fields_getpoint_f90(keyFlds_, iter.toFortran(), values.size(), values[0]);
  return oops::LocalIncrement(vars_, values, varlens);
}
// -----------------------------------------------------------------------------
void FieldsQG::setLocal(const oops::LocalIncrement & x, const GeometryQGIterator & iter) {
  const std::vector<double> vals = x.getVals();
  qg_fields_setpoint_f90(keyFlds_, iter.toFortran(), vals.size(), vals[0]);
}
// -----------------------------------------------------------------------------
size_t FieldsQG::serialSize() const {
  int nx, ny, nz, lbc;
  qg_fields_sizes_f90(keyFlds_, nx, ny, nz);
  qg_fields_lbc_f90(keyFlds_, lbc);
  size_t nn = nx * ny * nz;
  if (lbc == 1) {
    nn += 2 * (nx + 1) * nz;
  }
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void FieldsQG::serialize(std::vector<double> & vect)  const {
  size_t size_fld = this->serialSize() - time_.serialSize();

  // Allocate space for fld, xb and qb
  std::vector<double> v_fld(size_fld, 0);

  // Serialize the field
  qg_fields_serialize_f90(keyFlds_, static_cast<int>(size_fld), v_fld.data());
  vect.insert(vect.end(), v_fld.begin(), v_fld.end());

  // Serialize the date and time
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void FieldsQG::deserialize(const std::vector<double> & vect, size_t & index) {
  int indexInt = static_cast<int>(index);
  qg_fields_deserialize_f90(keyFlds_, static_cast<int>(vect.size()), vect.data(), indexInt);
  index = static_cast<size_t>(indexInt);
  time_.deserialize(vect, index);
}
// -----------------------------------------------------------------------------
}  // namespace qg
