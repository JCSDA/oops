/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsBiasIncrement.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/ObsBias.h"
#include "model/ObsBiasCovariance.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
ObsBiasIncrement::ObsBiasIncrement(const ObsSpaceQG &, const eckit::Configuration & conf)
  : bias_(ObsBias::ntypes, 0.0), active_(ObsBias::ntypes, false)
{
  const eckit::LocalConfiguration covconf = conf.getSubConfiguration("covariance");

  active_[0] = covconf.has("stream");
  active_[1] = covconf.has("uwind");
  active_[2] = covconf.has("vwind");
  active_[3] = covconf.has("wspeed");
  bool on = false;
  std::string strn = "";
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (jj > 0) strn += ", ";
    if (active_[jj]) {
      strn += "on";
      on = true;
    } else {
      strn += "off";
    }
  }
  if (on) {oops::Log::trace() << "ObsBiasIncrement created : " << strn << std::endl;}
}
// -----------------------------------------------------------------------------
ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const bool copy)
  : bias_(ObsBias::ntypes, 0.0), active_(other.active_)
{
  if (copy) {
    for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] = other.bias_[jj];
  }
  this->makePassive();
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::makePassive() {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (!active_[jj]) bias_[jj] = 0.0;
  }
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    bias_[jj] = b1[jj] - b2[jj];
  }
  this->makePassive();
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::zero() {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] = 0.0;
}
// -----------------------------------------------------------------------------
ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] = rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] += rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] -= rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] *= fact;
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] += fact * rhs.bias_[jj];
  this->makePassive();
}
// -----------------------------------------------------------------------------
double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (active_[jj]) zz += bias_[jj] * rhs.bias_[jj];
  }
  return zz;
}
// -----------------------------------------------------------------------------
double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  int ii = 0;
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (active_[jj]) {
      zz += bias_[jj] * bias_[jj];
      ++ii;
    }
  }
  if (ii > 0) zz = std::sqrt(zz/ii);
  return zz;
}
// -----------------------------------------------------------------------------
size_t ObsBiasIncrement::serialSize() const {
  size_t nn = bias_.size();
  return nn;
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::serialize(std::vector<double> & vect) const {
  vect.insert(vect.end(), bias_.begin(), bias_.end());
  oops::Log::trace() << "ObsBiasIncrement::serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::deserialize(const std::vector<double> & vect, size_t & index) {
  for (unsigned int jj = 0; jj < bias_.size(); ++jj) {
    bias_[jj] = vect[index];
    ++index;
  }
  oops::Log::trace() << "ObsBiasIncrement::deserialize done" << std::endl;
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::print(std::ostream & os) const {
  bool on = false;
  std::string strn = "";
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (jj > 0) strn += ", ";
    if (active_[jj]) {
      on = true;
      std::ostringstream strs;
      strs << bias_[jj];
      strn += strs.str();
    } else {
      strn += "0.0";
    }
  }
  if (on) os << std::endl << "ObsBiasIncrement = " << strn;
}
// -----------------------------------------------------------------------------
}  // namespace qg
