/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <math.h>

#include "oops/util/Logger.h"

#include "model/ObsDataQG.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgFortran.h"

#include "eckit/exception/Exceptions.h"

namespace qg {
// -----------------------------------------------------------------------------
ObsVecQG::ObsVecQG(const ObsSpaceQG & obsdb, const std::string & name)
  : obsdb_(obsdb), keyOvec_(0)
{
  qg_obsvec_setup_f90(keyOvec_, obsdb.assimvariables().size(), obsdb.nobs());
  if (!name.empty()) obsdb_.getdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVecQG::ObsVecQG(const ObsVecQG & other)
  : obsdb_(other.obsdb_), keyOvec_(0) {
  qg_obsvec_clone_f90(keyOvec_, other.keyOvec_);
  qg_obsvec_copy_f90(keyOvec_, other.keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVecQG::~ObsVecQG() {
  qg_obsvec_delete_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVecQG & ObsVecQG::operator= (const ObsVecQG & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  qg_obsvec_copy_f90(keyOvec_, keyOvecRhs);
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
ObsVecQG & ObsVecQG::operator=(const ObsDataQG<float> & rhs) {
  *this = rhs.vect();
  return *this;
}
// -----------------------------------------------------------------------------
void ObsVecQG::zero() {
  qg_obsvec_zero_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecQG::setToMissing(int ii) {
  qg_obsvec_settomissing_ith_f90(keyOvec_, ii);
}
// -----------------------------------------------------------------------------
void ObsVecQG::ones() {
  qg_obsvec_ones_f90(keyOvec_);
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
  qg_obsvec_random_f90(obsdb_, keyOvec_);
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
  int iobs;
  qg_obsvec_nobs_f90(keyOvec_, iobs);
  double zz = 0.0;
  if (iobs > 0) {
    qg_obsvec_dotprod_f90(keyOvec_, keyOvec_, zz);
    zz = sqrt(zz/iobs);
  }
  return zz;
}
// -----------------------------------------------------------------------------
void ObsVecQG::mask(const ObsVecQG & mask) {
  qg_obsvec_mask_with_missing_f90(keyOvec_, mask.toFortran());
}
// -----------------------------------------------------------------------------
void ObsVecQG::save(const std::string & name) const {
  obsdb_.putdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
Eigen::VectorXd ObsVecQG::packEigen(const ObsVecQG & mask) const {
  Eigen::VectorXd vec(packEigenSize(mask));
  qg_obsvec_get_withmask_f90(keyOvec_, mask.toFortran(), vec.data(), vec.size());
  return vec;
}
// -----------------------------------------------------------------------------
size_t ObsVecQG::packEigenSize(const ObsVecQG & mask) const {
  int nobs;
  qg_obsvec_nobs_withmask_f90(keyOvec_, mask.toFortran(), nobs);
  return nobs;
}
// -----------------------------------------------------------------------------
void ObsVecQG::read(const std::string & name) {
  obsdb_.getdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVecQG::print(std::ostream & os) const {
  if (nobs() == 0) {
    os << obsdb_.obsname() << " no observations.";
  } else {
    double zmin, zmax, zavg;
    qg_obsvec_stats_f90(keyOvec_, zmin, zmax, zavg);
    std::ios_base::fmtflags f(os.flags());
    os << obsdb_.obsname() << " nobs= " << nobs()
       << "  Min=" << zmin
       << ", Max=" << zmax
       << ", Average=" << zavg;
    os.flags(f);
  }
}
// -----------------------------------------------------------------------------
unsigned int ObsVecQG::nobs() const {
  int iobs;
  qg_obsvec_nobs_f90(keyOvec_, iobs);
  unsigned int nobs(iobs);
  return nobs;
}
// -----------------------------------------------------------------------------
size_t ObsVecQG::size() const {
  int iobs;
  qg_obsvec_size_f90(keyOvec_, iobs);
  size_t nobs(iobs);
  return nobs;
}
// -----------------------------------------------------------------------------
}  // namespace qg
