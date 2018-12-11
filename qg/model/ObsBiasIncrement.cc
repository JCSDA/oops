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
ObsBiasIncrement::ObsBiasIncrement(const eckit::Configuration & conf)
  : bias_(ObsBias::ntypes, 0.0), active_(ObsBias::ntypes, false)
{
  active_[0] = conf.has("stream");
  active_[1] = conf.has("uwind");
  active_[2] = conf.has("vwind");
  active_[3] = conf.has("wspeed");
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
ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const eckit::Configuration &)
  : bias_(ObsBias::ntypes, 0.0), active_(other.active_)
{
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] = other.bias_[jj];
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
void ObsBiasIncrement::serialize(std::vector<double> & vect) const {
  double s_obs_bias = bias_.size();
  double s_obs_act = active_.size();
  vect.push_back(s_obs_bias + s_obs_act + 2);
  vect.push_back(s_obs_bias);
  vect.insert(vect.end(), bias_.begin(), bias_.end());
  vect.push_back(s_obs_act);
  vect.insert(vect.end(), active_.begin(), active_.end());
  oops::Log::trace() << "ObsBiasIncrement::serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::deserialize(const std::vector<double> & vect) {
  unsigned int s_obs_bias = std::lround(vect[0]);
  for (unsigned int jj = 0; jj < s_obs_bias; ++jj) bias_[jj] = vect[jj + 1];
  unsigned int s_obs_act = std::lround(vect[s_obs_bias + 1]);
  for (unsigned int jj = 0; jj < s_obs_act; ++jj) active_[jj] = vect[jj + s_obs_bias + 2];
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
