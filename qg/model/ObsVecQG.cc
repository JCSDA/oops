/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <math.h>

#include "util/Logger.h"

#include "model/ObsVecQG.h"
#include "model/ObsSpaceQG.h"
#include "model/QgFortran.h"

namespace qg {
// -----------------------------------------------------------------------------
ObsVecQG::ObsVecQG(const ObsSpaceQG & obsdb)
  : obsdb_(obsdb), keyOvec_(0)
{
  qg_obsvec_setup_f90(keyOvec_, obsdb.nout(), obsdb.nobs());
}
// -----------------------------------------------------------------------------
ObsVecQG::ObsVecQG(const ObsVecQG & other, const bool copy)
  : obsdb_(other.obsdb_), keyOvec_(0) {
  qg_obsvec_clone_f90(other.keyOvec_, keyOvec_);
  if (copy) {
    qg_obsvec_assign_f90(keyOvec_, other.keyOvec_);
  } else {
    qg_obsvec_zero_f90(keyOvec_);
  }
}
// -----------------------------------------------------------------------------
ObsVecQG::~ObsVecQG() {
  qg_obsvec_delete_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVecQG & ObsVecQG::operator= (const ObsVecQG & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  qg_obsvec_assign_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecQG & ObsVecQG::operator*= (const double & zz) {
  qg_obsvec_mul_scal_f90(keyOvec_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecQG & ObsVecQG::operator+= (const ObsVecQG & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  qg_obsvec_add_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecQG & ObsVecQG::operator-= (const ObsVecQG & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  qg_obsvec_sub_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecQG & ObsVecQG::operator*= (const ObsVecQG & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  qg_obsvec_mul_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVecQG & ObsVecQG::operator/= (const ObsVecQG & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  qg_obsvec_div_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
void ObsVecQG::zero() {
  qg_obsvec_zero_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecQG::axpy(const double & zz, const ObsVecQG & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  qg_obsvec_axpy_f90(keyOvec_, zz, keyOvecRhs);
}
// -----------------------------------------------------------------------------
void ObsVecQG::invert() {
  qg_obsvec_invert_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecQG::random() {
  qg_obsvec_random_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
double ObsVecQG::dot_product_with(const ObsVecQG & other) const {
  const int keyOvecOther = other.keyOvec_;
  double zz;
  qg_obsvec_dotprod_f90(keyOvec_, keyOvecOther, zz);
  return zz;
}
// -----------------------------------------------------------------------------
double ObsVecQG::rms() const {
  double zz;
  qg_obsvec_dotprod_f90(keyOvec_, keyOvec_, zz);
  int iobs;
  qg_obsvec_nobs_f90(keyOvec_, iobs);
  zz = sqrt(zz/iobs);
  return zz;
}
// -----------------------------------------------------------------------------
void ObsVecQG::read(const std::string & name) {
  obsdb_.getdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecQG::save(const std::string & name) const {
  obsdb_.putdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecQG::print(std::ostream & os) const {
  double zmin, zmax, zavg;
  qg_obsvec_minmaxavg_f90(keyOvec_, zmin, zmax, zavg);
  os << obsdb_.obsname() << " nobs= " << size()
     << " Min=" << zmin << ", Max=" << zmax << ", Average=" << zavg;
}
// -----------------------------------------------------------------------------
unsigned int ObsVecQG::size() const {
  int iobs;
  qg_obsvec_nobs_f90(keyOvec_, iobs);
  unsigned int nobs(iobs);
  return nobs;
}
// -----------------------------------------------------------------------------
}  // namespace qg
