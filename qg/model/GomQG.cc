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
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace qg {

// -----------------------------------------------------------------------------
GomQG::GomQG(const LocationsQG & locs, const oops::Variables & vars,
             const std::vector<size_t> & sizes):
  vars_(vars)
{
  // gom_setup just creates and allocates the GeoVaLs object without filling
  // in values
  qg_gom_setup_f90(keyGom_, locs, vars_);
}
// -----------------------------------------------------------------------------
/*! QG GeoVaLs Constructor with Config */

  GomQG::GomQG(const eckit::Configuration & config,
               const ObsSpaceQG & ospace, const oops::Variables & vars):
  vars_(vars)
{
  qg_gom_create_f90(keyGom_, vars_);
  qg_gom_read_file_f90(keyGom_, config);
}
// -----------------------------------------------------------------------------
// Copy constructor
GomQG::GomQG(const GomQG & other):
  vars_(other.vars_)
{
  qg_gom_create_f90(keyGom_, vars_);
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
void GomQG::read(const eckit::Configuration & config) {
  qg_gom_read_file_f90(keyGom_, config);
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
  os << " nobs= " << nobs << std::scientific << std::setprecision(4)
     << "  Min=" << std::setw(12) << zmin
     << ", Max=" << std::setw(12) << zmax
     << ", RMS=" << std::setw(12) << zrms;
  os.flags(f);

  // If the min value across all variables is positive, then this may be an
  // error measurement.  If so, print the location and variable where the
  // maximum occurs to the debug stream, for use in debugging

  if (zmin >= 0.0) {
    double mxval;
    int iloc;
    oops::Variables maxvar;

    qg_gom_maxloc_f90(keyGom_, mxval, iloc, maxvar);

    oops::Log::debug() << "GomQG: Maximum Value = " << std::setprecision(4)
                       << mxval << " at location = " << iloc
                       << " and variable = " << maxvar << std::endl;
  }
}
// -----------------------------------------------------------------------------
}  // namespace qg
