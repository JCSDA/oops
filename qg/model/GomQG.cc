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

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "util/Logger.h"
#include "model/LocationsQG.h"
#include "model/QgFortran.h"
#include "model/VariablesQG.h"

namespace qg {

// -----------------------------------------------------------------------------
GomQG::GomQG(const LocationsQG & locs, const oops::Variables & var) {
  const VariablesQG varqg(var);
  qg_gom_setup_f90(keyGom_, locs.toFortran(), varqg.toFortran());
}
// -----------------------------------------------------------------------------
GomQG::GomQG(const eckit::Configuration & config) {
  qg_gom_create_f90(keyGom_);
  const eckit::Configuration * conf = &config;
  qg_gom_read_file_f90(keyGom_, &conf);
}
// -----------------------------------------------------------------------------
GomQG::~GomQG() {
  qg_gom_delete_f90(keyGom_);
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
GomQG & GomQG::operator*=(const double & zz) {
  qg_gom_mult_f90(keyGom_, zz);
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
  const eckit::Configuration * conf = &config;
  qg_gom_read_file_f90(keyGom_, &conf);
}
// -----------------------------------------------------------------------------
void GomQG::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  qg_gom_write_file_f90(keyGom_, &conf);
}
// -----------------------------------------------------------------------------
void GomQG::print(std::ostream & os) const {
  int nn;
  double zmin, zmax, zavg;
  qg_gom_minmaxavg_f90(keyGom_, nn, zmin, zmax, zavg);
  os << " nobs= " << nn << " Min=" << zmin << ", Max=" << zmax << ", RMS=" << zavg;
}
// -----------------------------------------------------------------------------
}  // namespace qg
