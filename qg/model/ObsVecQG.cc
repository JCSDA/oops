/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <math.h>
#include <algorithm>
#include <sstream>

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
void ObsVecQG::serialize(std::vector<double> & values) const {
  int nn;
  qg_obsvec_size_f90(keyOvec_, nn);
  std::vector<double> newvalues(nn);
  qg_obsvec_serialize_f90(keyOvec_, newvalues.data(), newvalues.size());
  values.reserve(values.size() + nn);
  values.insert(values.end(), newvalues.begin(), newvalues.end());
}
// -----------------------------------------------------------------------------
void ObsVecQG::deserialize(const std::vector<double> & values, size_t & ind) {
  int indInt = static_cast<int>(ind);
  qg_obsvec_deserialize_f90(keyOvec_, values.data(), static_cast<int>(values.size()), indInt);
  ind = static_cast<size_t>(indInt);
}
// -----------------------------------------------------------------------------
void ObsVecQG::maskAndSerialize(const ObsVecQG & mask, std::vector<double> & values) const {
  int nobs;
  qg_obsvec_nobs_withmask_f90(keyOvec_, mask.toFortran(), nobs);
  std::vector<double> newvalues(nobs);
  qg_obsvec_get_withmask_f90(keyOvec_, mask.toFortran(), newvalues.data(), newvalues.size());
  values.insert(values.end(), newvalues.begin(), newvalues.end());
}
// -----------------------------------------------------------------------------
std::vector<size_t> ObsVecQG::maskAndSerialIndices(const ObsVecQG & mask) const {
  int nobs;
  qg_obsvec_nobs_withmask_f90(keyOvec_, mask.toFortran(), nobs);
  std::vector<size_t> indices(nobs);
  std::vector<int> indices_int(nobs);
  qg_obsvec_getindices_withmask_f90(keyOvec_, mask.toFortran(),
                                    indices_int.data(), indices_int.size());
  std::transform(indices_int.begin(), indices_int.end(), indices.begin(),
                   [](int ele) { return static_cast<size_t>(ele); });
  return indices;
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
    os << std::left << std::setw(8) << obsdb_.obsname() << std::right
       << std::setw(0) << " nobs= " << std::setw(5) << nobs()
       << std::setw(0) << "  Min="  << std::setw(12) << zmin
       << std::setw(0) << ", Max="  << std::setw(12) << zmax
       << std::setw(0) << ", Average=" << std::setw(12) << zavg;
    os.flags(f);
  }
}
// -----------------------------------------------------------------------------
std::string ObsVecQG::info(const std::string & prefix) const {
  std::stringstream ss;
  this->print(ss);
  std::string grep = "\n" + prefix;
  if (!grep.empty() && std::isalnum(grep.back())) grep += ": ";
  return "\n" + ss.str();
}
// -----------------------------------------------------------------------------
std::string ObsVecQG::info(const std::string & prefix, const ObsDataQG<int> &) const {
  std::stringstream ss;
  this->print(ss);
  std::string grep = "\n" + prefix;
  if (!grep.empty() && std::isalnum(grep.back())) grep += ": ";
  return grep + ss.str();
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
size_t ObsVecQG::serialSize() const {
  return this->size();
}
// -----------------------------------------------------------------------------
void ObsVecQG::readAppended(const std::string & name) {
  throw eckit::NotImplemented("ObsVecQG::readAppended() is not implemented.", Here());
}
// -----------------------------------------------------------------------------
void ObsVecQG::zeroAppended() {
  throw eckit::NotImplemented("ObsVecQG::zeroAppended() is not implemented.", Here());
}
// -----------------------------------------------------------------------------
}  // namespace qg
