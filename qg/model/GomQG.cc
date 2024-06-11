/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/GomQG.h"

#include <iomanip>

#include "eckit/config/Configuration.h"
#include "model/LocationsQG.h"
#include "model/ObsSpaceQG.h"
#include "model/QgFortran.h"
#include "model/QgTraitsFwd.h"
#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace qg {

// -----------------------------------------------------------------------------
GomQG::GomQG(const Locations_ & locs, const oops::Variables & vars,
             const std::vector<size_t> & sizes):
  vars_(vars)
{
// The QG model cannot handle locations sampled with more than one method yet.
  ASSERT(locs.numSamplingMethods() == 1);

// All variables have same levels
  for (size_t jj = 1; jj < sizes.size(); ++jj) ASSERT(sizes[jj] == sizes[0]);
  const int levs = sizes[0];
  qg_gom_setup_f90(keyGom_, locs.samplingMethod(0).sampledLocations().size(), vars_, levs);
}
// -----------------------------------------------------------------------------
/*! QG GeoVaLs Constructor with Config */

  GomQG::GomQG(const eckit::Configuration & config,
               const ObsSpaceQG & ospace, const oops::Variables & vars):
  vars_(vars)
{
  qg_gom_create_f90(keyGom_);
  qg_gom_read_file_f90(keyGom_, vars_, config);
}
// -----------------------------------------------------------------------------
// Copy constructor
GomQG::GomQG(const GomQG & other):
  vars_(other.vars_)
{
  qg_gom_create_f90(keyGom_);
  qg_gom_copy_f90(keyGom_, other.keyGom_);
}
// -----------------------------------------------------------------------------
GomQG::~GomQG() {
  qg_gom_delete_f90(keyGom_);
}
// -----------------------------------------------------------------------------
double GomQG::rms() const {
  double zz;
  qg_gom_rms_f90(keyGom_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
double GomQG::normalizedrms(const GomQG & rhs) const {
  GomQG temp_GomQG(*this);
  qg_gom_divide_f90(temp_GomQG.keyGom_, rhs.keyGom_);
  return temp_GomQG.rms();
}
// -----------------------------------------------------------------------------
void GomQG::zero() {
  qg_gom_zero_f90(keyGom_);
}
// -----------------------------------------------------------------------------
void GomQG::random() {
  qg_gom_random_f90(keyGom_);
}
// -----------------------------------------------------------------------------
GomQG & GomQG::operator=(const GomQG & rhs) {
  const int keyGomRhs = rhs.keyGom_;
  qg_gom_copy_f90(keyGom_, keyGomRhs);
  return *this;
}
// -----------------------------------------------------------------------------
GomQG & GomQG::operator*=(const double & zz) {
  qg_gom_mult_f90(keyGom_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
GomQG & GomQG::operator+=(const GomQG & other) {
  qg_gom_add_f90(keyGom_, other.keyGom_);
  return *this;
}
// -----------------------------------------------------------------------------
GomQG & GomQG::operator-=(const GomQG & other) {
  qg_gom_diff_f90(keyGom_, other.keyGom_);
  return *this;
}
// -----------------------------------------------------------------------------
GomQG & GomQG::operator*=(const GomQG & other) {
  qg_gom_schurmult_f90(keyGom_, other.keyGom_);
  return *this;
}
// -----------------------------------------------------------------------------
double GomQG::dot_product_with(const GomQG & other) const {
  double zz;
  qg_gom_dotprod_f90(keyGom_, other.keyGom_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void GomQG::fill(const oops::Variable &var, const ConstVectorRef<size_t> &indx,
                 const ConstMatrixRef<double> &vals, const bool) {
  const size_t npts = indx.size();
  const size_t nlev = vals.cols();
  std::vector<int> findx(indx.size());
  for (Eigen::Index jj = 0; jj < indx.size(); ++jj) findx[jj] = indx[jj] + 1;

  qg_gom_fill_f90(keyGom_, var.name().size(), var.name().data(), npts,
                  findx.data(), nlev, vals.data());
}
// -----------------------------------------------------------------------------
void GomQG::fillAD(const oops::Variable &var, const ConstVectorRef<size_t> &indx,
                   MatrixRef<double> vals, const bool) const {
  const size_t npts = indx.size();
  const size_t nlev = vals.cols();
  std::vector<int> findx(indx.size());
  for (Eigen::Index jj = 0; jj < indx.size(); ++jj) findx[jj] = indx[jj] + 1;

  qg_gom_fillad_f90(keyGom_, var.name().size(), var.name().data(), npts,
                    findx.data(), nlev, vals.data());
}
// -----------------------------------------------------------------------------
void GomQG::read(const eckit::Configuration & config) {
  qg_gom_read_file_f90(keyGom_, vars_, config);
}
// -----------------------------------------------------------------------------
void GomQG::write(const eckit::Configuration & config) const {
  qg_gom_write_file_f90(keyGom_, config);
}
// -----------------------------------------------------------------------------
void GomQG::print(std::ostream & os) const {
  int nobs;
  double zmin, zmax, zrms;
  qg_gom_stats_f90(keyGom_, nobs, zmin, zmax, zrms);
  std::ios_base::fmtflags f(os.flags());
  os << " nobs= " << nobs
     << "  Min=" << zmin
     << ", Max=" << zmax
     << ", RMS=" << zrms;
  os.flags(f);
}
// -----------------------------------------------------------------------------
}  // namespace qg
