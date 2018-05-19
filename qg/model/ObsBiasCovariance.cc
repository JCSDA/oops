/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsBiasCovariance.h"

#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "util/Logger.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
ObsBiasCovariance::ObsBiasCovariance(const eckit::Configuration & conf)
  : conf_(conf), variance_(ObsBias::ntypes, 0.0)
{
//  if (!conf.empty()) {
    std::vector<double> zz(4, 0.0);
    if (conf.has("stream")) zz[0] = conf.getDouble("stream");
    if (conf.has("uwind"))  zz[1] = conf.getDouble("uwind");
    if (conf.has("vwind"))  zz[2] = conf.getDouble("vwind");
    if (conf.has("wspeed")) zz[3] = conf.getDouble("wspeed");

    std::string strn = "";
    for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
      if (jj > 0) strn += ", ";
      if (std::abs(zz[jj]) > 1.0e-8) {
        variance_[jj] = zz[jj] * zz[jj];
        std::ostringstream strs;
        strs << variance_[jj];
        strn += strs.str();
      } else {
        variance_[jj] = 0.0;
        strn += "0.0";
      }
    }
    Log::info() << "ObsBiasCovariance created, variances = " << strn << std::endl;
//  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::multiply(const ObsBiasIncrement & dxin,
                                 ObsBiasIncrement & dxout) const {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dxout[jj] = dxin[jj] * variance_[jj];
    } else {
      dxout[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::inverseMultiply(const ObsBiasIncrement & dxin,
                                        ObsBiasIncrement & dxout) const {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dxout[jj] = dxin[jj] / variance_[jj];
    } else {
      dxout[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::randomize(ObsBiasIncrement & dx) const {
  static std::mt19937 generator(4);
  static std::normal_distribution<double> distribution(0.0, 1.0);
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      double zz = distribution(generator);
      dx[jj] = zz * std::sqrt(variance_[jj]);
    } else {
      dx[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::print(std::ostream & os) const {
  os << "ObsBiasCovariance::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace qg
